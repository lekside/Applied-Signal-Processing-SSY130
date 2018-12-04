# Simple script to convert soundfont and midi files to c-compatible arrays

from os import listdir
from textwrap import fill


def file2c(fn, array_name):
    """ Read a file and return a C-compatible array with its binary data and its length """

    # Open file and generate a list with hex-encoded byte contents
    with open(fn, 'rb') as f:
        rawdata = f.read()
    hexdata = [hex(x) for x in rawdata]

    # Generate C-compatible array

    # Header
    c_str = 'const char ' + array_name + '[' + str(len(hexdata)) + '] = {'

    # Data contents
    for x in hexdata:
        c_str += x + ', '

    # Remove last ', ' elements
    c_str = c_str[:-2]

    # Tail
    c_str += '};'
    return c_str, len(hexdata)


def files2c(readdir, ext):
    """ Read all files in a given directory with a given extension and convert them to C-compatible arrays """
    names = listdir(readdir)

    # Generate list of files matching target extension
    files = []
    for name in names:
        if ext == name.split('.')[-1]:
            files.append(name)

    # Generate all C-arrays
    c_strs = []
    lens = []
    for file in files:
        str, len = file2c(readdir + '/' + file, 'mp3_data')
        c_strs.append(str)
        lens.append(len)

    return c_strs, files, lens

# Convert all .mp3 files in an ./raw directory to binary C-arrays
# Write output to mp3_data.h/.c

# Read relevant data from files

datas, files, lens = files2c('./raw', 'mp3')

# Generate output data

files_s = [file.split('.')[0] for file in files]


h_data =  '/** @file Contains raw mp3 data. Automatically generated by converter.py */\n'
h_data += '#ifndef MP3_DATA_H_\n'
h_data += '#define MP3_DATA_H_\n'
h_data += '#include "config.h"\n'
for idx, file in enumerate(files_s):
    h_data += '#ifdef MP3_DATA_{}\n'.format(file.upper())
    h_data += 'extern const char mp3_data[{}];\n'.format(lens[idx])
    h_data += '#endif\n'
h_data += '#endif /* MP3_DATA_{}_H_ */\n'.format(files_s[idx])

h_name = 'mp3_data.h'.format(files_s[idx])
c_name = 'mp3_data.c'.format(files_s[idx])

c_hdr = '#include "mp3_data.h"\n'

with open(h_name, 'w') as f:
    f.write(h_data)

with open(c_name, 'w') as f:
    f.write(c_hdr)
    for idx, file in enumerate(files_s):
        f.write('#ifdef MP3_DATA_{}\n'.format(file.upper()))
        f.write(fill(datas[idx], width=120, break_long_words=False, break_on_hyphens=False) + '\n')
        f.write('#endif\n')

