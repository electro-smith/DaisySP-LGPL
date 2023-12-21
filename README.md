<div align=center>
<img width=15% src="https://raw.githubusercontent.com/electro-smith/daisysp/master/resources/assets/banner.png">

# DaisySP-LGPL • LGPL DSP Modules for [DaisySP](https://www.github.com/electro-smith/DaisySP/)

[![Style Badge](https://github.com/electro-smith/DaisySP-LGPL/workflows/Style/badge.svg)](https://github.com/electro-smith/DaisySP-LPGL/actions?query=workflow%3AStyle)
[![Documentation Badge](https://github.com/electro-smith/DaisySP-LGPL/workflows/Documentation/badge.svg)](https://electro-smith.github.io/DaisySP-LGPL/index.html)
![Discord](https://img.shields.io/discord/1037767234803740694?logo=discord&label=Discord)
[![Forum Badge](https://img.shields.io/badge/chat-daisy%20forum-orange)](https://forum.electro-smith.com/)
[![License Badge](https://img.shields.io/badge/license-LGPL-yellow)](https://opensource.org/licenses/LGPL)

> DaisySP-LGPL is a set of LGPL licensed DSP (Digital Signal Processing) modules for use in the DaisySP library.

</div>
<br>

## 🧐 DaisySP and Contribution

This is not a standalone library, and requires [DaisySP](https://www.github.com/electro-smith/DaisySP) to work properly.

Any contributions to this library are welcome under the same guidelines as DaisySP.

## ⚠️ License

DaisySP-LGPL uses the LGPL-2.1 license.

For the full license, read the [LICENSE](https://github.com/electro-smith/DaisySP-LGPL/blob/master/LICENSE) file in the root directory.

## ✉️ Distribution
If you distribute code that makes use of DaisySP-LGPL, you are obligated to give end users the option to recompile with the LGPL components having been replaced.

To make this easier, we have provided a script called [gather_lgpl.sh](https://github.com/electro-smith/DaisySP-LGPL/tree/master/distribution/gather_lgpl.sh) which can gather the required files with a simple build system for that purpose.

1. From your project directory run `<DAISYSP-LGPL_DIR>/distribution/gather_lgpl.sh <PROJECT_NAME>`
2. Use the `-l` and `-d` flags to override the default libDaisy and DaisySP locations.
3. A standalone `distribution` folder will be created containing:
    - A Readme
    - The LICENSE file(s)
    - A standalone Makefile
    - And any relevant build files under a `resources` folder.
4. End users with the `distribution` folder can relink from source with their own modified copy of `libdaisy-lgpl.a` by replacing that file in the `resources` folder, and running `make`.
