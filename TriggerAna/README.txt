# e1039-analysis/TriggerAna
A package to study e1039 trigger
Copied and modified Kenichi's instructive module "e1039-analysis/AnaSimDst"


## Compilation

The analysis code is located in `e1039-analysis/TriggerAna/src`.
The following commands set up the shell environment and compile the code;
```
cd e1039-analysis/TriggerAna
source setup.sh
mkdir build
cd build
cmake-trigger-ana
make install
```

## Analysis Flow
0. Run simulation using e1039-analysis/SimChainDev/Fun4Sim.C
   - It generates several root files including DST.root. We need this DST file.
   - A little modified codes that I used are under "simulation"

1. DST -> Micro DST (uDST)
   - To reduce the data size (uDST has reduced number of branches).
   - The following command lists up all DSTs into `list_dst.txt`;
   ```
   cd e1039-analysis/TriggerAna/work
   ./make-dst-list.sh /your/dst/dir

   ex)
   ./make-dst-list.sh /pnfs/e906/persistent/users/kimmj/SimChainDev/E1039/*/out

   ```

   - The following commands create uDST from listed DSTs.
   - The uDST file is `uDST.root`.
   ```
   cd e1039-analysis/TriggerAna
   source setup.sh
   cd work
   root -b -q Fun4SimDst.C
   ```


2. uDST -> Simple Tree
   - To create a simple tree.
   - No event selection is made yet.	
   - To further reduce the data size & calculate your custom quantities.
   - The TTree object is outputted to `sim_tree.root`.
   ```
   root -b -q Fun4SimMicroDst.C
   ```

3. Simple Tree -> Road Tree
   - SpinQuest values are converted to road ids, etc.
   - For each event, road IDs of mu+ and mu- are saved. road ID=0 for out of acceptance muon
   ```
   root -b -q run_GetRoad.C

   	  


## Author

Minjung Kim <kminjung@umich.edu>

