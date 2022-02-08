## Simple Instructions for Testing 3D AMR Lossy Compressor

### Step 1 Download SZ

Download and install [SZ](https://github.com/szcompressor/SZ) into **SZDir** followed by the instructions provided by SZ's repo.

### Step 2: Download Data

Download and unzip data or use example data

**Method 1:** Download all datasets using [gdown](https://pypi.org/project/gdown/)
- gdown can be easily installed by pip or pip3 install, e.g., pip install gdown.
- Run download-data.sh bash sript.
```
cd data
./download-data.sh
```

**Method 2:** Download separate datasets directly from Google Drive
- Run1_Z10: https://drive.google.com/file/d/118ataZqVlCFq5D04MhtTTOBBvVO6DR6B/view?usp=sharing
- Run1_Z5: https://drive.google.com/file/d/1PHMUV6c9K5V9AN-TF_CYf9ROR86oZYSm/view?usp=sharing
- Run1_Z3: https://drive.google.com/file/d/1hEjnw3llHj2P9z3WyrWSQAWqPAQ-LKxy/view?usp=sharing
- Run1_Z2: https://drive.google.com/file/d/1xWwf2_YFDjod0mUIRl1iffj95ZemR4Q6/view?usp=sharing
- Run2_t2: https://drive.google.com/file/d/1xZnts_QRAN5Av3NbBzLzt0FybVafV7JC/view?usp=sharing
- Run2_t3: https://drive.google.com/file/d/1vy3Bry1RmAdXqMNXhwXFWZOzF6bCa7hM/view?usp=sharing
- Rune_t4: https://drive.google.com/file/d/1CeI_nxhU5_rrT2e78Jsy025t5lOmaQAd/view?usp=sharing

**Step 3.1: Build Our Compressor**
```
cd src
g++ amrcompressor.cpp -o amrcompressor -I SZDir/include/ -L SZDir/lib/ -lSZ -lzstd -lzlib -O3
```

**Step 3.2: Build Baseline & zMesh**
```
cd baseline
g++ naive1D.cpp -o naive1D -I SZDir/include/ -L SZDir/lib/ -lSZ -lzstd -lzlib -O3
g++ 3dBaseline.cpp -o 3dBaseline -I SZDir/include/ -L SZDir/lib/ -lSZ -lzstd -lzlib -O3
g++ zMesh.cpp -o zMesh -I SZDir/include/ -L SZDir/lib/ -lSZ -lzstd -lzlib -O3
```

**Step 4.1: Test Our Compressor & Baseline**
```
./amrcompressor ../data/run2_t3.bin SZDir/example/sz.config eb_level_0 eb_level_1 eb_level_2 eb_level_3
```

For example, all AMR levels with the same absolute error bound of 5E+9 for the example Nyx data:
```
./amrcompressor ../data/run2_t3.bin SZDir/example/sz.config 5E+9 5E+9 5E+9 5E+9
```

**Step 4.2: Test zMesh**

zMesh can only be run on **Run1**, will be modified for Run2 later.
```
./zMesh ../data/grid_z10.bin order256.bin SZDir/example/sz.config eb_level_0 eb_level_1 eb_level_2 eb_level_3
```

For example, all AMR levels with the same absolute error bound of 5E+9 for the Run1_Z10:
```
./zMesh ../data/grid_z10.bin order256.bin SZDir/example/sz.config 5E+9 5E+9 5E+9 5E+9
```
