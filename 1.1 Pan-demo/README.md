# Pan-genome demo
## This section is a toy example for testing various tools and methods
### Design
![fig](https://github.com/WengangXbio/OR_project/blob/14daadff4f5e3260c48add3ef91ee0fe9ff97277/1.1%20Pan-demo/pic/pan-design.png)
Three genome sequences (G1,G2,G3) are generated. Each segment is 200bp in length.
### Sequence design 
#### 1. Simulating SV
G1(S4) / G2(S5)
#### 2. Simulating deletion
G1(S2) / G2(*)
#### 3. Simulating SNP mutation
G1(S3)/ G2(S3*)
#### 4. Simulating insertion mutation (3bp)
G1(S1)/ G3(S1*)
#### 5. Simulating deletion mutation (3bp)
G1(S7)/ G3(S7*)

### 1. Results
#### 1. minigraph
![fig](https://github.com/WengangXbio/OR_project/blob/6fc20b0d83caee740d8a53b2c69e64d7f939217a/1.1%20Pan-demo/pic/minigraph.demo.png)
#### 2. PGGB
![fig](https://github.com/WengangXbio/OR_project/blob/6fc20b0d83caee740d8a53b2c69e64d7f939217a/1.1%20Pan-demo/pic/pggb.demo.png)

#### minigraph cannot produce correct graph, therefore PGGB is selected for following studies as trails

### 2. odgi functions
