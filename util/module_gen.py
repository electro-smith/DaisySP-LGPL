import datetime
import os
import re

"""
Generate module templates for the Daisy library.
This can be used as a CLI tool or a module.
Functions:
    template(moduleName, author) -> (headerString, cppString)
    createFiles(moduleName, author) -> files
"""

header = r"""
#pragma once
#ifndef DSY_{}_H
#define DSY_{}_H

namespace daicsp
{{

enum state
{{
    {}_STATE_A,
    {}_STATE_B,
    {}_LAST,
}};

/** {}
 *  Author: {}
 *  Date: {}
 */
class {}
{{
    public:
        {} () {{}}
        ~{} () {{}}

        /** Initializes the {} module.
         *  \param p - description
         */
        void Init();

        /** Processes a single sample and returns it.
         *  \param in - input sample
         */
        float Process(const float &in);

        /** Setter for private member param
         */
        inline void SetParam(const float &param) {{ param_ = param; }}

        /** Setter for private member complex_param
         */
        void SetComplexParam(const float &complex_param);

    private:
        float param_;
        float foo_bar_;
        float a_, b_;
}};

}} // namespace daicsp

#endif // DSY_{}_H
"""

cpp = """
#include "{}.h"

using namespace daicsp;

void {}::Init()
{{
    // Set private members to defaults
    param_ = 0.0f;
    a_ = 0.0f;
    b_ = 1.0f;

    // Do stuff
    
}}

float {}::Process(const float &in)
{{
    // Do something and return the output.
    return (in * param_) + a_ - b_;
}}

void {}::SetComplexParam(const float &complex_param)
{{
    a_ = complex_param;
    b_ = 1.0f - complex_param;
}}
"""

def _queryUser(prompt):
    response = input(prompt + ' [y/N] ')
    return response.lower().strip() in ('y', 'yes')

disk = re.compile(r'[A-Za-z]:[/\\]?')

# necessary because os.makedirs gets confused by '..'
def _makeDirs(path):
    # TODO -- make this compatible with windows separators (i.e. \)
    folders = path.split('/')
    tempPath = ''
    init = 0
    
    if folders[0] == '':
        tempPath += '/'
        init = 1
    elif disk.search(folders[0]) is not None:
        tempPath += folders[0]
        init = 1
    
    for folder in folders[init:]:
        tempPath = os.path.join(tempPath, folder)
        if not os.path.exists(tempPath):
            os.mkdir(tempPath)
    



def template(moduleName, author):
    """
    Generate template strings for the header and cpp files. The date is automatically
    generated from the system clock.
    """
    modUpper = moduleName.upper()
    modLower = moduleName.lower()
    headerData = [
        modUpper, modUpper, modUpper, modUpper, modUpper,
        moduleName, author, str(datetime.date.today()),
        moduleName, moduleName, moduleName, moduleName,
        modUpper,
    ]
    cppData = [
        modLower, moduleName, moduleName, moduleName
    ]
    return (header.lstrip('\n').format(*headerData), cpp.lstrip('\n').format(*cppData))


def createFiles(moduleName, author, path='.'):
    """
    Generate the header and cpp files. If the files already exist, the user
    is asked if they would like to overwrite them.
    """
    strings = template(moduleName, author)
    modLower = moduleName.lower()

    headerName = os.path.join(path, modLower + '.h')
    cppName = os.path.join(path, modLower + '.cpp')

    if not os.path.exists(path):
        if not _queryUser(f'The path {path} does not exist. Create?'):
            print('No files written.')
            return
        _makeDirs(path)
        
    if os.path.exists(headerName) or os.path.exists(cppName):
        if not _queryUser('One or both files already exist. Overwrite existing files?'):
            print('No files written.')
            return

    with open(headerName, 'w') as file:
        file.write(strings[0])
    with open(cppName, 'w') as file:
        file.write(strings[1])

    print('Successfully generated files for {}.'.format(moduleName))


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='Generate module templates for the Daisy library.'
    )
    parser.add_argument(
        'name', metavar='Module', type=str,
        help='Name of the module to be created, given in CamelCase.'
    )
    parser.add_argument(
        '-a', metavar='author', type=str, default='author',
        help='Author of the module.'
    )
    parser.add_argument(
        '-d', metavar='directory', type=str, default='.',
        help='Path to parent directory.'
    )

    args = parser.parse_args()
    createFiles(args.name, args.a, path=args.d)