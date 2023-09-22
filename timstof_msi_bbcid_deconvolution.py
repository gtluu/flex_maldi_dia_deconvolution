import argparse
import os
import sys
import subprocess
import numpy as np
import pandas as pd
import deimos
from pyteomics import mgf as pyt
from pyTDFSDK import *
from pyTDFSDK.classes import TdfData


def get_args():
    """
    Parse command line parameters.

    :return: Arguments with default or user specified values.
    :rtype: dict
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--ms1',
                        help='Path to MALDI-TIMS-MS .d directory containing analysis.tdf.',
                        required=True,
                        type=str)
    parser.add_argument('--ms2',
                        help='Path to MALDI-TIMS-MS/MS in bbCID mode .d directory containing analysis.tdf.',
                        required=True,
                        type=str)
    parser.add_argument('--outdir',
                        help='Path to directory where the output MGF file should be saved.',
                        required=True,
                        type=str)
    parser.add_argument('--outfile',
                        help='Filename of output MGF file.',
                        required=True,
                        type=str)
    parser.add_argument('--intensity_threshold',
                        help='Intensity threshold used when to subset spectrum from a frame. Defaults to 100. Set to 0 '
                             'to disable subset.',
                        default=100,
                        type=int)
    parser.add_argument('--corr_cutoff',
                        help='Correlation score cutoff between 0 and 1 to pass to DEIMoS for deconvolution workflow. '
                             'Defaults to 0.9.',
                        default=0.9,
                        type=float)
    parser.add_argument('--mz_low_tol',
                        help='Low end of the mass tolerance used during DEIMoS deconvolution workflow. Defaults to '
                             '-0.05 Da.',
                        default=-0.05,
                        type=float)
    parser.add_argument('--mz_high_tol',
                        help='High end of the mass tolerance used during DEIMoS deconvolution workflow. Defaults to '
                             '0.05 Da.',
                        default=0.05,
                        type=float)
    parser.add_argument('--mob_low_tol',
                        help='Low end of the mobility tolerance used during DEIMoS deconvolution workflow. Defaults to '
                             '-0.05 Da.',
                        default=-0.05,
                        type=float)
    parser.add_argument('--mob_high_tol',
                        help='High end of the mobility tolerance used during DEIMoS deconvolution workflow. Defaults '
                             'to 0.05 Da.',
                        default=0.05,
                        type=float)
    parser.add_argument('--mz_relative',
                        help='Flag to indicate mz_low_tol and mz_high_tol are in ppm instead of Daltons.',
                        action='store_true')
    parser.add_argument('--mob_relative',
                        help='Flag to indicate mob_low_tol and mob_high_tol are in percentage instead of 1/K0.',
                        action='store_true')
    parser.add_argument('--mz_res',
                        help='Resolution in the mass dimension. Defaults to 0.05 Da.',
                        default=0.05,
                        type=float)
    parser.add_argument('--mob_res',
                        help='Resolution in the mobility dimension. Defaults to 0.05 1/K0',
                        default=0.05,
                        type=float)
    parser.add_argument('--bin_size',
                        help='Bin size for binning deconvoluted MS/MS spectra before writing to MGF. Defaults to 0.1 '
                             'Da.',
                        default=0.1,
                        type=float)
    parser.add_argument('--min_peaks',
                        help='Minimum number of peaks required for deconvoluted MS/MS spectrum to be written to MGF '
                             'file. Defaults to 4.',
                        default=4,
                        type=int)
    return vars(parser.parse_args())


# Copied from TIMSCONVERT parse submodule.
def extract_3d_tdf_spectrum(tdf_data, frame, scan_begin, scan_end):
    """
    Extract spectrum from TDF data with m/z and intensity arrays. Spectrum can either be centroid or quasi-profile
    mode. "Raw" and "centroid" modes uses pyTDFSDK.tims.tims_read_scans_v2(). "Profile" mode data is not available due
    to the resulting data size.

    :param tdf_data: tdf_data object containing metadata from analysis.tdf database.
    :type tdf_data: timsconvert.classes.TimsconvertTdfData
    :param frame: Frame ID from the Frames table in analysis.tdf/analysis.tsf database.
    :type frame: int
    :param scan_begin: Beginning scan number (corresponding to 1/K0 value) within frame.
    :type scan_begin: int
    :param scan_end: Ending scan number (corresponding to 1/K0 value) within frame (non-inclusive).
    :type scan_end: int
    :return: Tuple of mz_array (np.array), intensity_array (np.array), and mobility_array (np.array) or
        (None, None, None) if spectra are empty.
    :rtype: tuple[numpy.array | None]
    """
    list_of_scans = tims_read_scans_v2(tdf_data.api, tdf_data.handle, frame, scan_begin, scan_end)
    frame_mz_arrays = []
    frame_intensity_arrays = []
    frame_mobility_arrays = []
    if scan_begin != 0:
        scan_end = scan_end - scan_begin
        scan_begin = 0
    for scan_num in range(scan_begin, scan_end):
        if list_of_scans[scan_num][0].size != 0 \
                and list_of_scans[scan_num][1].size != 0 \
                and list_of_scans[scan_num][0].size == list_of_scans[scan_num][1].size:
            mz_array = tims_index_to_mz(tdf_data.api, tdf_data.handle, frame, list_of_scans[scan_num][0])
            intensity_array = list_of_scans[scan_num][1]
            mobility = tims_scannum_to_oneoverk0(tdf_data.api, tdf_data.handle, frame, np.array([scan_num]))[0]
            mobility_array = np.repeat(mobility, mz_array.size)
            frame_mz_arrays.append(mz_array)
            frame_intensity_arrays.append(intensity_array)
            frame_mobility_arrays.append(mobility_array)
    if frame_mz_arrays and frame_intensity_arrays and frame_mobility_arrays:
        frames_array = np.stack((np.concatenate(frame_mz_arrays, axis=None),
                                 np.concatenate(frame_intensity_arrays, axis=None),
                                 np.concatenate(frame_mobility_arrays, axis=None)),
                                axis=-1)
        frames_array = np.unique(frames_array[np.argsort(frames_array[:, 0])], axis=0)
        mz_array = frames_array[:, 0]
        intensity_array = frames_array[:, 1]
        mobility_array = frames_array[:, 2]
        return mz_array, intensity_array, mobility_array
    else:
        return None, None, None


def load_data(args):
    """
    Load data from .d directory containing analysis.tdf SQLite database using pyTDFSDK.classes.TdfData.

    :return: Tuple of MS1 and MS/MS dataset.
    :rtype: tuple[pyTDFSDK.classes.TdfData]
    """
    bruker_dll = init_tdf_sdk_api()

    print('Reading MS1 Data')
    ms1 = TdfData(args['ms1'], bruker_dll)
    print('Reading MS/MS (bbCID) Data')
    ms2 = TdfData(args['ms2'], bruker_dll)

    if ms1.analysis['Frames'].size != ms2.analysis['Frames'].size:
        print('MS1 and MS/MS data should have the exact same number of frames and dimensions.')
        sys.exit(1)
    else:
        return ms1, ms2


def load_frame(ms, frame, intensity_threshold=100):
    """
    Load spectrum for a given frame using function copied from timsconvert.parse into a pd.DataFrame.

    :param ms: Dataset to pull spectrum from.
    :type ms:  pyTDFSDK.classes.TdfData
    :param frame: ID of the frame from Frames table to pull spectrum from.
    :type frame: int
    :param intensity_threshold: Intensity threshold used when to subset spectrum from a frame. Defaults to 100. Set to
        0 to disable subset.
    :return: Spectrum with m/z, intensity, and mobility values.
    :rtype: pd.DataFrame
    """
    ms_frames_dict = ms.analysis['Frames'][ms.analysis['Frames']['Id'] == frame].to_dict(orient='records')[0]
    ms_mz, ms_intensity, ms_mobility = extract_3d_tdf_spectrum(ms, frame, 0, int(ms_frames_dict['NumScans']))
    df = pd.DataFrame(data={'mz': ms_mz,
                            'inverse_reduced_ion_mobility': ms_mobility,
                            'intensity': ms_intensity})
    # Recommended to remove peaks with an absolute intensity of less than 100 counts.
    # Done to improve deconvolution time at the cost of losing low abundance peaks.
    # Without this step, the script is very RAM hungry.
    if intensity_threshold != 0:
        df = df[df['intensity'] >= intensity_threshold].reset_index()
    return df


def deconvolute_frame(ms1_df, ms2_df, corr_cutoff=0.9,
                      mz_low_tol=-0.05, mz_high_tol=0.05, mz_relative=False, mz_res=0.05,
                      mob_low_tol=-0.05, mob_high_tol=0.05, mob_relative=False, mob_res=0.05):
    """
    Deconvolute MS/MS spectrum for a given frame using DEIMoS.

    :param ms1_df: MS1 pd.DataFrame from load_frame().
    :type ms1_df:  pd.DataFrame
    :param ms2_df: MS/MS pd.DataFrame from load_frame().
    :type ms2_df:  pd.DataFrame
    :param corr_cutoff: Correlation score cutoff between 0 and 1 to pass to DEIMoS for deconvolution workflow. Defaults
        to 0.9.
    :type corr_cutoff: float
    :param mz_low_tol: Low end of the mass tolerance used during DEIMoS deconvolution workflow. Defaults to -0.05 Da.
    :type mz_low_tol: float
    :param mz_high_tol: High end of the mass tolerance used during DEIMoS deconvolution workflow. Defaults to 0.05 Da.
    :type mz_high_tol: float
    :param mz_relative: Indicate mz_low_tol and mz_high_tol are in ppm instead of Daltons. Defaults to False.
    :type mz_relative: bool
    :param mz_res: Resolution in the mass dimension. Defaults to 0.05 Da.
    :type mz_res: float
    :param mob_low_tol: Low end of the mobility tolerance used during DEIMoS deconvolution workflow. Defaults to -0.05
        Da.
    :type mob_low_tol: float
    :param mob_high_tol: High end of the mobility tolerance used during DEIMoS deconvolution workflow. Defaults to 0.05
        Da.
    :type mob_high_tol: float
    :param mob_relative: Indicate mob_low_tol and mob_high_tol are in percentage instead of 1/K0. Defaults to False.
    :type mob_relative: bool
    :param mob_res: Resolution in the mobility dimension. Defaults to 0.05 1/K0.
    :type mob_res: float
    :return: pandas object grouped by precursor m/z values.
    :rtype: pandas.core.groupby.DataFrameGroupBy
    """
    res = deimos.deconvolution.deconvolve_ms2(ms1_df, ms1_df, ms2_df, ms2_df,
                                              cluster_kwargs={'dims': ['inverse_reduced_ion_mobility'],
                                                              'tol': [mob_high_tol],
                                                              'relative': [mob_relative]},
                                              profile_kwargs={
                                                  'dims': ['mz', 'inverse_reduced_ion_mobility'],
                                                  'low': [mz_low_tol, mob_low_tol],
                                                  'high': [mz_high_tol, mob_high_tol],
                                                  'relative': [mz_relative, mob_relative],
                                                  'resolution': [mz_res, mob_res]},
                                              apply_kwargs={'dims': ['inverse_reduced_ion_mobility']})
    res = res.loc[(res['inverse_reduced_ion_mobility_score'] > corr_cutoff)]
    res = res.groupby(by=[x for x in res.columns if x.endswith('_ms1')], as_index=False).agg(list)
    return res


def write_spectrum_to_mgf(args, row, count, low_mz, high_mz, bin_size=0.1, min_peaks=4):
    """
    Write deconvoluted MS/MS spectrum to MGF file using pyteomics.mgf.

    :param args: Arguments with default or user specified values.
    :type args: dict
    :param row: Row from a pd.DataFrame containing deconvoluted spectra from DEIMoS as a dict.
    :type row: dict
    :param count: Scan count used for scan number in MGF file.
    :type count: int
    :param low_mz: Low m/z range used for spectrum binning.
    :type low_mz: float
    :param high_mz: High m/z range used for spectrum binning.
    :type high_mz: float
    :param bin_size: Bin size for binning deconvoluted MS/MS spectra before writing to MGF. Defaults to 0.1 Da.
    :type bin_size: float
    :param min_peaks: Minimum number of peaks required for deconvoluted MS/MS spectrum to be written to MGF file.
        Defaults to 4.
    """
    mz_array = np.array(row['mz_ms2'], dtype=np.float64)
    intensity_array = np.array(row['intensity_ms2'], dtype=np.float32)
    sorted_index = mz_array.argsort()
    mz_array = mz_array[sorted_index]
    intensity_array = intensity_array[sorted_index]

    # Round m/z to 4 decimal places and group m/z values by summing intensities.
    bin_df = pd.DataFrame({'mz': mz_array, 'intensity': intensity_array})
    bin_df = bin_df.round({'mz': 4})
    bin_df = bin_df.groupby(by='mz', as_index=False).aggregate(sum)
    mz_array = bin_df['mz']
    intensity_array = bin_df['intensity']

    # Bin deconvoluted MS/MS spectra before writing to MGF.
    bins = np.arange(low_mz, high_mz, bin_size, dtype=np.float64)
    unique_indices, inverse_indices = np.unique(np.digitize(mz_array, bins), return_inverse=True)
    bin_counts = np.bincount(inverse_indices)
    np.place(bin_counts, bin_counts < 1, [1])
    mz_array = np.bincount(inverse_indices, weights=mz_array) / bin_counts
    intensity_array = np.bincount(inverse_indices, weights=intensity_array)

    # Only write MS/MS spectrum if there are N or more fragment peaks.
    if mz_array.size >= min_peaks:
        ms2_dict = {'m/z array': mz_array,
                    'intensity array': intensity_array,
                    'params': {'FEATURE_ID': count,
                               'PEPMASS': float(row['mz_ms1']),
                               'SCANS': count,
                               'MSLEVEL': 2,
                               'ION_MOBILITY': float(row['inverse_reduced_ion_mobility_ms1']),
                               'CHARGE': '1+'}}
        pyt.write([ms2_dict], output=os.path.join(args['outdir'], args['outfile']) + '.mgf', file_mode='a')


def main():
    """
    Run workflow.
    """
    args = get_args()
    print(args)
    ms1, ms2 = load_data(args)

    count = 1
    for frame in range(1, ms1.analysis['Frames'].shape[0]):
        print('Parsing Frame: ' + str(frame))
        df1 = load_frame(ms1, frame, args['intensity_threshold'])
        df2 = load_frame(ms2, frame, args['intensity_threshold'])

        if not df1.empty and not df2.empty:
            try:
                print('Deconvoluting Frame: ' + str(frame))
                res = deconvolute_frame(df1, df2, args['corr_cutoff'],
                                        args['mz_low_tol'], args['mz_high_tol'],
                                        args['mz_relative'], args['mz_res'],
                                        args['mob_low_tol'], args['mob_high_tol'],
                                        args['mob_relative'], args['mob_res'])

                print('Writing MGF for Coordinates Frame ' + str(frame))
                for index, row in res.iterrows():
                    write_spectrum_to_mgf(args,
                                          row,
                                          count,
                                          int(float(ms1.analysis['GlobalMetadata']['MzAcqRangeLower'])),
                                          int(float(ms1.analysis['GlobalMetadata']['MzAcqRangeUpper'])),
                                          args['bin_size'],
                                          args['min_peaks'])
                    count += 1
            except ValueError:
                # Skip frame if ValueError occurs, most likely due to low quality spectra for that frame.
                print('Value Error')
                pass


if __name__ == '__main__':
    main()
