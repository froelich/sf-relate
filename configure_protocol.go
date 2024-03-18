package relativeMatch

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"time"

	"github.com/BurntSushi/toml"
	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/ldsec/lattigo/v2/ckks"
	"go.dedis.ch/onet/v3/log"
)

type BasicProtocolConfig struct {
	NumMainParties int `toml:"num_main_parties"`
	HubPartyId     int `toml:"hub_party_id"`

	CkksParams string `toml:"ckks_params"`

	divSqrtMaxLen   int `toml:"div_sqrt_max_len"`
	Servers         map[string]mpc.Server
	MpcFieldSize    int `toml:"mpc_field_size"`
	MpcDataBits     int `toml:"mpc_data_bits"`
	MpcFracBits     int `toml:"mpc_frac_bits"`
	MpcNumThreads   int `toml:"mpc_num_threads"`
	LocalNumThreads int `toml:"local_num_threads"`

	SharedKeysPath string `toml:"shared_keys_path"`
	BindingIP      string `toml:"binding_ipaddr"`
	Para           int    `toml:"PARA"`
}

type BasicProtocolInfo struct {
	MpcObj mpc.ParallelMPC
	Cps    *crypto.CryptoParams
	Config *BasicProtocolConfig
}

func InitializeBasicProtocol(pid int, configFolder string, network mpc.ParallelNetworks) (relativeProt *BasicProtocolInfo) {

	config := new(BasicProtocolConfig)
	// Import global parameters
	if _, err := toml.DecodeFile(filepath.Join(configFolder, "configGlobal.toml"), config); err != nil {
		fmt.Println(err)
		return
	}

	var chosen int
	switch config.CkksParams {
	case "PN12QP109":
		chosen = ckks.PN12QP109
	case "PN13QP218":
		chosen = ckks.PN13QP218
	case "PN14QP438":
		chosen = ckks.PN14QP438
	case "PN15QP880":
		chosen = ckks.PN15QP880
	case "PN16QP1761":
		chosen = ckks.PN16QP1761
	default:
		panic("Undefined value of CKKS params in Config")
	}

	params := ckks.DefaultParams[chosen]
	prec := uint(config.MpcFieldSize)
	var networks mpc.ParallelNetworks
	if network == nil {
		log.LLvl1(config.BindingIP)
		log.LLvl1(config.Servers)
		log.LLvl1(config.NumMainParties)
		log.LLvl1(config.MpcNumThreads)
		networks = mpc.ParallelNetworks(mpc.InitCommunication(config.BindingIP, config.Servers, pid, config.NumMainParties+1, config.MpcNumThreads, config.SharedKeysPath))
	} else {
		networks = network
	}

	for thread := range networks {
		networks[thread].SetMHEParams(params)
	}
	var rtype mpc_core.RElem
	switch config.MpcFieldSize {
	case 256:
		rtype = mpc_core.LElem256Zero
	case 128:
		rtype = mpc_core.LElem128Zero
	default:
		panic("Unsupported value of MPC field size")
	}

	log.LLvl1(fmt.Sprintf("MPC parameters: bit length %d, data bits %d, frac bits %d",
		config.MpcFieldSize, config.MpcDataBits, config.MpcFracBits))
	mpcEnv := mpc.InitParallelMPCEnv(networks, rtype, config.MpcDataBits, config.MpcFracBits)
	for thread := range mpcEnv {
		mpcEnv[thread].SetHubPid(config.HubPartyId)
	}

	cps := networks.CollectiveInit(params, prec)
	// Access shared random seed to use for sampling hashing parameters
	// The seed is sampled from the Network[0]
	networks[0].Rand.SwitchPRG(-1)
	seed := make([]byte, 16)
	networks[0].Rand.RandRead(seed)
	networks[0].Rand.RestorePRG()
	log.LLvl1("Seed for shared hashing randomness:", seed)

	// save the seed in a file named "seed.bin" to be read by the python scrip0t
	file, err := os.OpenFile(filepath.Join(config.SharedKeysPath, "seed.bin"), os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		panic(err)
	}
	defer file.Close()
	// write the bytes in seed to the file
	_, err = file.Write(seed)
	if err != nil {
		panic(err)
	}

	if pid != 0 {
		// Run MHE step 0 to generate the shared hashing randomness
		log.LLvl1("Running step 0 to generate shared parameters")
		out, err := exec.Command("python3", "notebooks/step0_sample_shared_keys.py", "-PARTY", strconv.Itoa(pid), "-FOLDER", configFolder).Output()
		if err != nil {
			log.LLvl1(string(out))
			panic(err)
		}
		log.LLvl1("Finished sampling shared hash table parameters!")
		out, err = exec.Command("python3", "notebooks/step1_hashing.py", "-PARTY", strconv.Itoa(pid), "-FOLDER", configFolder).Output()
		log.LLvl1("Done hashing.")
		if err != nil {
			log.LLvl1(string(out))
			panic(err)
		}
		log.LLvl1(string(out))
	}
	// need to synchronize
	// wait for all parties to finish
	synchronize(mpcEnv)

	log.LLvl1("Moving onto the actual MHE (step 2)")
	return &BasicProtocolInfo{
		MpcObj: mpcEnv,
		Cps:    cps,
		Config: config,
	}
}

func synchronize(mpcEnv []*mpc.MPC) {
	mainMPCObj := mpcEnv[0]
	pid := mainMPCObj.GetPid()

	var dummy mpc_core.RElem = mainMPCObj.GetRType().Zero()
	if pid == 0 {
		for p := 1; p < mainMPCObj.GetNParty(); p++ {
			dummy = mainMPCObj.Network.ReceiveRElem(dummy, p)
		}
		for p := 1; p < mainMPCObj.GetNParty(); p++ {
			mainMPCObj.Network.SendRData(dummy, p)
		}
	} else {
		mainMPCObj.Network.SendRData(dummy, 0)
		_ = mainMPCObj.Network.ReceiveRElem(dummy, 0)
	}
}

func (g *ProtocolInfo) SyncAndTerminate(closeChannelFlag bool) {
	synchronize(g.basicProt.MpcObj)

	if closeChannelFlag {
		// Close all threads
		for t := range g.basicProt.MpcObj {
			g.basicProt.MpcObj[t].Network.CloseAll()
		}
		// need to sleep for a bit to make sure all channels are closed on all parties, before doing anything else
		time.Sleep(2 * time.Second)
	}
}
