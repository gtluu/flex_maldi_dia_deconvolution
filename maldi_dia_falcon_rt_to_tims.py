import os
import argparse


def get_args():
    """
    Parse command line parameters.

    :return: Arguments with default or user specified values.
    :rtype: dict
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',
                        help='Path to directory where the output MGF file should be saved.',
                        required=True,
                        type=str)
    return vars(parser.parse_args())


input = get_args()['input']

# Change RTINSECONDS tag back to ION_MOBILITY.
with open(input, 'r') as mgf_file:
    mgf = mgf_file.read()
mgf = mgf.replace('RTINSECONDS', 'ION_MOBILITY')
with open(os.path.join(os.path.dirname(input),
                       os.path.splitext(os.path.split(input)[1])[0])[:-4] + '.mgf',
          'w') as mgf_file:
    mgf_file.write(mgf)

os.remove(input[:-18] + '_tmp.mgf')
os.remove(input)
