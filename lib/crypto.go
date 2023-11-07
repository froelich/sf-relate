package lib

import (
	"math"

	"github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/ldsec/lattigo/v2/ckks"
	"go.dedis.ch/onet/v3/log"
)

func CSubPlain(cryptoParams *crypto.CryptoParams, X crypto.PlainVector, Y crypto.CipherVector) crypto.CipherVector {
	res := make(crypto.CipherVector, len(X))
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range Y {
			res[i] = eval.SubNew(X[i], Y[i])
		}
		return nil
	})
	return res

}

func CSub(cryptoParams *crypto.CryptoParams, X crypto.CipherVector, Y crypto.CipherVector) crypto.CipherVector {
	res := make(crypto.CipherVector, len(X))
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range Y {
			res[i] = eval.SubNew(X[i], Y[i])
		}
		return nil
	})
	return res

}

// changed this Other thing to actually implement Plain in Y and Cipher in X
func CSubPlainOther(cryptoParams *crypto.CryptoParams, X crypto.CipherVector, Y crypto.PlainVector) crypto.CipherVector {
	res := make(crypto.CipherVector, len(X))
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range X {
			res[i] = eval.SubNew(X[i], Y[i])
		}
		return nil
	})
	return res

}

func CRescale(cryptoParams *crypto.CryptoParams, X crypto.CipherVector) crypto.CipherVector {
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range X {
			eval.Rescale(X[i], cryptoParams.Params.Scale(), X[i])
		}
		return nil
	})
	return X

}

func InvSqrt(x complex128) complex128 {
	return complex(1.0/math.Sqrt(real(x)), 0)
}

func TestFunc(x complex128) complex128 {
	return complex(real(x)*real(x)*99+real(x)*real(x)*real(x), 0)
}

func CInvSqrtSquareApprox(sourcePid int, net *mpc.Network, cryptoParams *crypto.CryptoParams, ctIn crypto.CipherVector, intv crypto.IntervalApprox, netDec *mpc.Network) crypto.CipherVector {
	res := make(crypto.CipherVector, len(ctIn))
	for i := range ctIn {
		res[i] = InvSqrtSquareApprox(sourcePid, net, cryptoParams, ctIn[i], intv, netDec)
	}
	return res
}

// InvSqrtApprox computes an encrypted approximated version of the inverse function
func InvSqrtSquareApprox(sourcePid int, net *mpc.Network, cryptoParams *crypto.CryptoParams, ctIn *ckks.Ciphertext,
	intv crypto.IntervalApprox, netDec *mpc.Network) *ckks.Ciphertext {
	var y *ckks.Ciphertext

	if intv.Degree == 0 {
		// this code is hanging locally for some reason
		ctDecrypt := net.CollectiveDecrypt(cryptoParams, ctIn, sourcePid)
		cdDecode := crypto.DecodeFloatVector(cryptoParams, crypto.PlainVector{ctDecrypt})
		cdOut := make([]float64, len(cdDecode))
		log.LLvl1("before approx result ", cdDecode[:300])
		for i := 0; i < len(cdDecode); i++ {
			cdOut[i] = 1.0 / math.Sqrt(cdDecode[i])
		}
		log.LLvl1("approx result ", cdOut[:300])
		cdOutEncrypted, _ := crypto.EncryptFloatVector(cryptoParams, cdOut)
		return cdOutEncrypted[0]
	}
	if !intv.InverseNew {
		cheby := ckks.Approximate(InvSqrt, complex(intv.A, 0), complex(intv.B, 0), intv.Degree)
		// this demonstrate the order of the coefficients are correct
		// We evaluate the interpolated Chebyshev interpolant on y
		// Change of variable
		a := cheby.A()
		b := cheby.B()

		if ctIn.Level() < int(math.Ceil(math.Log2(float64(intv.Degree+1)))+1)+3 {
			// log.LLvl1("scale bef boot ", ctIn.Level(), ctIn.Scale())
			net.CollectiveBootstrap(cryptoParams, ctIn, sourcePid)
			// log.LLvl1("scale aft boot ", ctIn.Level(), ctIn.Scale())
		}

		// seems like the tranformation from the comment was incorrect, this is correct
		err := cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			y = eval.MultByConstNew(ctIn.CopyNew().Ciphertext(), 2/(b-a))
			if err := eval.Rescale(y, cryptoParams.Params.Scale(), y); err != nil {
				panic(err)
			}
			eval.AddConst(y, (-a-b)/(b-a), y)
			return nil
		})
		if err != nil {
			log.Fatal(err)
		}

		ctInCopy := y.CopyNew().Ciphertext()
		err = cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			var err error
			y, err = eval.EvaluateCheby(ctInCopy, cheby, ctInCopy.Scale())
			if err != nil {
				log.Fatal(err)
			}
			return err
		})
		if err != nil {
			log.Fatal(err)
		}

		if y.Level() < cryptoParams.BootLimit+3 {
			// log.LLvl1("bootstrap before iter")
			net.CollectiveBootstrap(cryptoParams, y, sourcePid)

		}

		c0 := float64(1.875)
		c1 := float64(-1.250)
		c2 := float64(0.375)

		//ctIn = crypto.LevelTest(crypto.CipherVector{ctIn}, cryptoParams, 2, " ", "ctIn in invSqrt")[0]
		var xc1, x2 *ckks.Ciphertext
		err = cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			xc1 = eval.MultByConstNew(ctIn, c1)
			eval.Rescale(xc1, cryptoParams.Params.Scale(), xc1)
			x2 = eval.MulRelinNew(ctIn, ctIn)
			eval.Rescale(x2, cryptoParams.Params.Scale(), x2)
			return nil
		})
		if err != nil {
			log.Fatal(err)
		}

		for i := 0; i < intv.Iter; i++ {
			if y.Level() < cryptoParams.BootLimit+3 {
				net.CollectiveBootstrap(cryptoParams, y, sourcePid)
				// log.LLvl1("bootstrap beginning of loop, iter ", i, " level ", y.Level())
			}

			var y2, y4 *ckks.Ciphertext
			err = cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
				y2 = eval.MulRelinNew(y, y)
				err = eval.Rescale(y2, cryptoParams.Params.Scale(), y2)
				if err != nil {
					log.Fatal(err)
				}
				y4 = eval.MulRelinNew(y2, y2)
				err = eval.Rescale(y4, cryptoParams.Params.Scale(), y4)
				if err != nil {
					log.Fatal(err)
				}
				return nil
			})
			if err != nil {
				log.Fatal(err)
			}
			// log.LLvl1("After y2 y4: ", y.Level())

			err = cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
				tmp0 := eval.MulRelinNew(xc1, y2)
				err = eval.Rescale(tmp0, cryptoParams.Params.Scale(), tmp0)
				if err != nil {
					log.Fatal(err)
				}
				eval.AddConst(tmp0, c0, tmp0)
				eval.MulRelin(tmp0, y, tmp0)
				err = eval.Rescale(tmp0, cryptoParams.Params.Scale(), tmp0)
				if err != nil {
					log.Fatal(err)
				}

				tmp1 := eval.MultByConstNew(y, c2)
				err = eval.Rescale(tmp1, cryptoParams.Params.Scale(), tmp1)
				if err != nil {
					log.Fatal(err)
				}
				eval.MulRelin(tmp1, x2, tmp1)
				err = eval.Rescale(tmp1, cryptoParams.Params.Scale(), tmp1)
				if err != nil {
					log.Fatal(err)
				}
				eval.MulRelin(tmp1, y4, tmp1)
				err = eval.Rescale(tmp1, cryptoParams.Params.Scale(), tmp1)
				if err != nil {
					log.Fatal(err)
				}

				y = eval.AddNew(tmp0, tmp1)
				// log.LLvl1("After tmp: ", y.Level())
				return nil
			})
			if err != nil {
				log.Fatal(err)
			}
		}
	} else {
		log.Fatal("Not implemented")
	}

	return y
}

func FindMaxIndx(pt3 []complex128, imagine bool) (float64, int) {
	ptInd := 0
	ptMax := 0.0
	for i := 0; i < len(pt3); i++ {
		if imagine {
			if math.Abs(imag(pt3[i])) > ptMax {
				ptMax = math.Abs(imag(pt3[i]))
				ptInd = i
			}
		} else if math.Abs(real(pt3[i])) > ptMax {
			ptMax = math.Abs(real(pt3[i]))
			ptInd = i
		}
	}
	return ptMax, ptInd
}

func FindMaxIndxWithUp(pt3 []complex128, imagine bool, up float64) (float64, int) {
	ptInd := 0
	ptMax := 0.0
	for i := 0; i < len(pt3); i++ {
		if imagine {
			if math.Abs(imag(pt3[i])) > ptMax && math.Abs(imag(pt3[i])) < up {
				ptMax = math.Abs(imag(pt3[i]))
				ptInd = i
			}
		} else if math.Abs(real(pt3[i])) > ptMax && math.Abs(real(pt3[i])) < up {
			ptMax = math.Abs(real(pt3[i]))
			ptInd = i
		}
	}
	return ptMax, ptInd
}

func InnerSumPartial(cryptoParams *crypto.CryptoParams, X crypto.CipherVector, Xsize int) crypto.CipherVector {
	// Sum all the ciphertexts in vector (i.e vector sum)
	res := make(crypto.CipherVector, len(X))
	for i := range X {
		res[i] = X[i].CopyNew().Ciphertext()
		rt := X[i].CopyNew().Ciphertext()
		for rotate := 1; rotate < Xsize; rotate++ {
			log.LLvl1("Rotate by ", 1)
			cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
				rt = eval.RotateNew(rt, 1)
				eval.Add(res[i], rt, res[i])
				return nil
			})
		}
	}
	return res
}

func ConstRotAndAdd(cryptoParams *crypto.CryptoParams, X crypto.CipherVector, Xsize int, nbrofRepeatperCipher int) crypto.CipherVector {
	if Xsize > cryptoParams.Params.Slots() {
		cipherSize := Xsize / cryptoParams.Params.Slots()
		res := make(crypto.CipherVector, len(X)/cipherSize)
		for i := 0; i < len(X); i++ {
			res[i%cipherSize] = crypto.Add(cryptoParams, res[i%cipherSize], X[i])
		}
		return res
	} else {
		// Potential improvment: when Y takes less than a ciphertext
		//nbrOfVecPerCipher := cryptoParams.Params.Slots() / Xsize
		resTmp := crypto.CopyEncryptedVector(X)
		for i := 0; i < len(X); i++ {
			rt := X[i].CopyNew().Ciphertext()
			for rotate := 1; rotate < nbrofRepeatperCipher; rotate++ {
				log.LLvl1("other rotate by ", Xsize)
				cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
					rt = eval.RotateNew(rt, Xsize)
					eval.Add(rt, resTmp[i], resTmp[i])
					return nil
				})
			}
		}
		res := make(crypto.CipherVector, 1)
		for i := 0; i < len(X); i++ {
			if i == 0 {
				res[0] = resTmp[i]
			} else {
				res[0] = crypto.Add(cryptoParams, res[0], resTmp[i])
			}

		}
		return res
	}
}

func InitEncryptedVector(cryptoParams *crypto.CryptoParams, dy int, scale float64) (crypto.CipherVector, int) {
	vector := make([]float64, dy)

	return crypto.EncryptFloatVectorWithScale(cryptoParams, vector, scale)
}
