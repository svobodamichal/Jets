1) first step is to run embedding - simulated jets (SP/PYTHIA) are embedded into real event and delta-pT is calculated
submitter/submit_embedding.sh

2) from the delta-pT histograms the response matrices are created
submitter/submit_buildRM.sh

3) detector response matrices are produced using Toymodel
~/jet_analysis/toymodel/submitter/submit_buildRM.sh

4) full response matrix is created by multiplying the BG and detector matrices
response_matrix/multiply_matrix.sh

5) for each prior function a response matrix compatible with roounfold is created
submitter/submit_buildRMROO.sh BGD

