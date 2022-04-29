# LaCyTools
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/84d61117770f4f5782a1decfed54b9cb)](https://www.codacy.com/gh/Tarskin/LaCyTools/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Tarskin/LaCyTools&amp;utm_campaign=Badge_Grade)

LaCyTools is open-source software for researcher focussing on chromatgraphy coupled to mass spectrometry. It includes modules for retention time and m/z calibration, peak integration, quality control and more. The tool was first described in a 2016 Journal of Proteome Research article, located at https://pubs.acs.org/doi/10.1021/acs.jproteome.6b00171

* __Source__: https://github.com/Tarskin/LaCyTools
* __Bug reports__: https://github.com/Tarskin/LaCyTools/issues

LaCyTools requires at least Python 3.8 and depends on NumPy, SciPy and matplotlib. NumPy and SciPy are used predominantly for data fitting and numerical operations on arrays. Matplotlib is used to draw all of the graphical elements that are used/created by LaCyToolsTools in both the GUI and the images (see the requirements.txt for the full list)

# Developer information
If you would like to take part in LaCyTools development, take a look at the CONTRIBUTING and CODE_OF_CONDUCT files.

# License information
See the file LICENSE.txt for information on the history of this software, terms & conditions for usage, and a DISCLAIMER OF ALL WARRANTIES.

# Installing LaCyTools (Windows)

* Install Python3.8
* Download the LaCyTools zip file and extract it to a non-protected folder
* Create a Python3.8 virtual environment (python3.8 -m venv LaCyTools) in the root folder where LaCyTools is extracted
* Activate the virtual environment (source LaCyTools/bin/activate)
* Install the specific requirements (pip install -r requirements.txt)