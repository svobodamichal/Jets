# How To: AuAu Jet Analysis

## Repository Location
The code can be found at:

```
/gpfs01/star/pwg/svomich/Jets
```

It consists of two main components:
- **AssociatedTriggerJet**: Runs over data.
- **EmbeddingPico**: For embedding analysis.

Both components have similar structures.

---

## Setup Instructions

1. **Required Files to Copy** (for both data and embedding):
   - `StRoot`
   - `starSubmit`
   - `run14AuAu200GeVPrescales`
   - `BadRunList_14.list`
   - `pico_low_14.list`
   - `runLocal.csh`
   - `ULTIMATErun.sh`
   - `test.list`
   - `fastjet1`

2. **FastJet Setup**:
   Create a symbolic link:

   ```bash
   ln -s /gpfs01/star/pwg/svomich/Jets/AssociatedTriggerJets/fastjet1/fastjet_install/include/fastjet /gpfs01/star/pwg/svomich/Jets/AssociatedTriggerJets/fastjet
   ```
    - You need to do this also for the Embedding part

3. **FastJet Path in Code**:
   - Currently, `StPicoHFJetMaker` requires a hardcoded path for FastJet.
   - Update the FastJet address manually due to RCF issues with dynamic paths.

4. **Environment**:
   - Use `starver pro`
   - With the above components, the code should run smoothly.

---

## Running the Code

- **Local Testing (on a few picos)**:
  ```bash
  csh runLocal.csh
  ```
  > Adjust pT hard bin inside the script if needed.

- **Job Submission**:
  1. Create folders `jobs` and `production`.
  2. Run:
     ```bash
     ./ULTIMATErun.sh
     ```

- **Scripts for Further Analysis**:
  Located in:
  ```
  /Jets/scripts/ForSasha
  ```
---

## After Job Completion

1. **Merging Outputs**:
   - **Embedding**:
     ```bash
     ./movefiles.sh
     ```
   - **Data**:
     ```bash
     ./merge.sh <folder_name>
     ```

2. **Local Post-Processing**:
   - Download results for all pThat bins.
   - Collect files from `MatchEff`, `RM` folders.

3. **Combine and Rebin**:
   ```bash
   bash plotall.sh
   ```

4. **Matching Efficiency**:
   ```bash
   root ploteffi.C
   ```

5. **Prepare Data for Unfolding**:
   - Download merged data file.
   - Run:
     ```bash
     root preparehistogramsforunfolding.C
     ```
   - Located in `DataPreparation` folder.

---

## Unfolding Procedure

1. **Required Folder Structure in Main Directory**:
   ```
   matchingeffi/              # contains pythia6_effi_<suffix>.root
   raw_spectra_to_unfold/HT2/ # contains histos_HT2_trig_<suffix>.root
   responseM_HT2/Pythia6/     # contains pythia6_all_<suffix>.root
   ```

2. **Copy Required Code**:
   - Copy `StRoot` from main directory.
   - Code located at:
     ```
     /Jets/StRoot/macros/unfolding_testPythia6
     ```

3. **Set Base Path**:
   Edit the file:
   ```
   Jets/StRoot/macros/set_paths.sh
   ```
   - Change `BASE` to your working directory.

4. **Suffix Convention**:
   - Used to differentiate between analysis versions.
   - Set in `1_run_unfolding.sh`:
     ```bash
     SYS="Main"
     ```
   - All input files must match the suffix, e.g.:
     - `pythia6_effi_Main.root`
     - `histos_HT2_trig_Main.root`
     - `pythia6_all_Main.root`

5. **Run Unfolding**:
   ```bash
   ./1_run_unfolding.sh
   ```

6. **Post-Unfolding Plots**:
   - Output directory: `Unfolded_HT2_<suffix>`
   - Go to `/unfolding` and run:
     ```bash
     bash plotratios.sh
     ```
---
