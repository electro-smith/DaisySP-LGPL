<h1>
  <img width=3% src="https://raw.githubusercontent.com/electro-smith/daisysp/master/resources/assets/banner.png">
  DaisySP-LGPL • LGPL Licensed DSP Modules for <a href="https://www.github.com/electro-smith/DaisySP/">DaisySP</a>
</h1>

[![Build Badge](https://github.com/electro-smith/DaisySP-LGPL/workflows/Build/badge.svg)](https://github.com/electro-smith/DaisySP-LGPL/actions?query=workflow%3ABuild)
[![Style Badge](https://github.com/electro-smith/DaisySP-LGPL/workflows/Style/badge.svg)](https://github.com/electro-smith/DaisySP-LPGL/actions?query=workflow%3AStyle)
[![Documentation Badge](https://github.com/electro-smith/DaisySP-LGPL/workflows/Documentation/badge.svg)](https://electro-smith.github.io/DaisySP-LGPL/index.html)
![Discord](https://img.shields.io/discord/1037767234803740694?logo=discord&label=Discord)
[![Forum Badge](https://img.shields.io/badge/chat-daisy%20forum-orange)](https://forum.electro-smith.com/)
[![License Badge](https://img.shields.io/badge/license-LGPL-yellow)](https://opensource.org/licenses/LGPL)

> DaisySP-LGPL is a set of LGPL licensed DSP (Digital Signal Processing) modules for use in the DaisySP library.

## DaisySP and Contribution

This is not a standalone library, and requires [DaisySP](https://www.github.com/electro-smith/DaisySP) to work properly.

Any contributions to this library are welcome under the same guidelines as DaisySP.

## License

DaisySP uses the MIT license.

It can be used in both closed source and commercial projects, and does not provide a warranty of any kind.

For the full license, read the [LICENSE](https://github.com/electro-smith/DaisySP/blob/master/LICENSE) file in the root directory.

## Distribution
If you distribute code that makes use of DaisySP-LGPL, you are obligated to give end users the option to recompile with the LGPL components having been replaced.

To make this easier, we have provided a script called [gather_lgpl.sh](https://github.com/electro-smith/libDaisy/tree/master/core/gather_lgpl.sh) which can gather the required files with a simple build system for that purpose.

1. From your project directory run `<LIBDAISY_DIR>/core/gather_lgpl.sh <PROJECT_NAME>`
2. Use the `-l` and `-d` flags to override the default libDaisy and DaisySP locations.
3. A standalone `lgpl` folder will be created containing:
    - A Readme
    - A standalone Makefile
    - And any relevant build files under a `resources` folder.
4. End users with the `lgpl` folder can relink from source with their own modified copy of `libdaisy-lgpl.a` by replacing that file in the `resources` folder, and running `make`