```
gridrun_data.sh
gridsub_data.sh
```
These are scripts to run real data in grid with the macro, `RecoE1039Data.C`. Since `/data2/` area is not acessible during job submission, these scripts assume you copied the dst data files from `/data2/` area to `/pnfs/e906/persistent/users/$USER/DstRun` area. Please refer to the [wiki page](https://github.com/E1039-Collaboration/e1039-wiki/wiki/data-file-on-grid) about setting the file location. To run list of the runs from `run_list.txt` in grid, use the following command;
```
./gridsub_data.sh ouput_dir 1 run_list.txt 0 single 
```
Here, the "0" options lets you run all the available events in the run. For the longer runs which can't be finished within the 24 hours limit, you need to split the data in some spill number interval. The `splitDST.C`macro
splits the longer DST runs in specified interval of spills. For now, the default number of spill interval for splitting is set to 750. 

```
root -b -q splitDst.C\(run_number\)
```

The splitted dst files will be outputed to `/pnfs/e906/persistent/users/$USER/DstRun/run_number/` (You might see some error msg regarding the persmission to write in persistent area, please ignore the message for now). After splitting the data for a particular run, use the following command to run all the splitted files from the runs in run_list.txt;

```
./gridsub_data.sh ouput_dir 1 run_list.txt 0 split 
```

## Note:
1. Please turn on `reco->set_enable_eval_dst(true);` in the `RecoE1039Data.C` macro to add the tracklet vector node to DST
1. Don't the chage the name of the input file in the macro while running grid job. Keep it as it is (i.e. `data.root`)
