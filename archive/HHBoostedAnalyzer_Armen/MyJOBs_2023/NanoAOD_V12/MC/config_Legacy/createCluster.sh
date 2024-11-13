#!/bin/bash

flavour="longlunch"

if [ ! -d "output" ]; then
    mkdir output
fi
if [ ! -d "error" ]; then
    mkdir error
fi
if [ ! -d "log" ]; then
    mkdir log
fi

echo 'executable                = run.csh' > Run.sub
echo 'arguments                 = RootFile_$(ProcId)' >> Run.sub
echo 'should_transfer_files     = YES' >>Run.sub
echo 'transfer_input_files      = RootFile_$(ProcId)' >> Run.sub
echo 'output                    = output/run.$(ClusterId).$(ProcId).out' >> Run.sub
echo 'error                     = error/run.$(ClusterId).$(ProcId).err' >> Run.sub
echo 'log                       = log/run.$(ClusterId).log' >> Run.sub
echo '+JobFlavour               =' '"'$flavour'"' >> Run.sub
