# MALDI DIA MSI Deconvolution from timsTOF fleX Data

## Installation and Usage

1. Create a new conda venv running Python 3.8.
```
conda create -n msi_dia_workflow python=3.8
```
2. Activate the new conda venv.
```
conda activate msi_dia_workflow
```
3. Install dependencies for this workflow
```
pip install -r requirements.txt
```
4. Run ```timstof_msi_bbcid_deconvolution.py```. The workflow requires a path to an MS1 imaging dataset and bbCID MS/MS imaging dataset. Use the ```--help``` parameter to see required parameters.
```
python timstof_msi_bbcid_deconvolution.py --ms1 ms1_dataset.d --ms2 ms2_dataset.d --outdir ~/output_directory --outfile output_filename
```
5. Prepare MGF file for clustering by running ```maldi_dia_falcon_tims_to_rt.py```. This replaces the temporarily ```ION_MOBILITY``` field as a ```RTINSECONDS``` field for ```falcon-ms``` compatibility.
```
python maldi_dia_falcon_tims_to_rt.py --input ~/output_directory/output_filename.mgf
```
6. Create a second new conda venv running Python 3.8.
```
conda create -n falcon python=3.8
```
7. Activate and install dependencies for falcon-ms.
```
conda activate falcon
pip install -r falcon_requirements.txt
```
8. Run falcon-ms clustering. Parameters can be changed as needed. Use ```falcon-ms --help``` for a full list of parameters.
```
falcon ~/output_directory/output_filename_tmp.mgf ~/output_directory/output_filename_clustered_tmp --export_representatives --eps 0.7 --precursor_tol 0.1 Da --rt_tol 0.1 --fragment_tol 0.1 --min_mz 100.0 --max_mz 2000.0 --min_peaks 2 --overwrite
```
9. Run ```maldi_dia_falton_rt_tims.py``` to revert ```RTINSECONDS``` to ```ION_MOBILITY``` and clean up tmp files.
```
python maldi_dia_falton_rt_tims.py --input ~/output_directory/output_filename_clustered_tmp.mgf
```
