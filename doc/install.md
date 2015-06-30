# Installation of STRTprep

To run STRTprep, you need (i) OSX or Linux, (ii) at least 8GB RAM, and (iii) enough storage space - for example, 100GB for processing of two 48-plex STRT libraries and three lanes of HiSeq sequencing per library.

## Installation of prerequisite softwares

First of all, you need to install several softwares to use Homebrew/Linuxbrew package management framework before installation of STRTprep; STRTprep installs additional softwares via the Homebrew/Linuxbrew during the installation of STRTPrep.

### OSX

- Check [Requirements](https://github.com/Homebrew/homebrew/blob/master/share/doc/homebrew/Installation.md#requirements) of [Homebrew](https://github.com/Homebrew/homebrew), and install Command Line Tool for Xcode.

### Linux

- Check [Dependencies](https://github.com/Homebrew/linuxbrew#dependencies) of [Linuxbrew](https://github.com/Homebrew/linuxbrew).
- Installation of the following softwares might help you and Linuxbrew.
  - mysql-devel (for kent-tools on CentOS)
  - imake (for openssl on CentOS)
  - python-setuptools (for cmake on CentOS)
