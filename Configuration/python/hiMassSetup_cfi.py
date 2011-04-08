### This files controls whether we are working with data/MC, signalMC (for skim), and channel (not implemented yet).
### Make sure to set it to the desired process.

### signal refers to signal MC => data must be false when signal is true
### data refers to data/mc analysis
### channels: emu, etau, mutau, tautau
global signal
global data
global channel

signal = False
data = True
channel = "tautau"
