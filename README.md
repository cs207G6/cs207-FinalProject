[![Build Status](https://travis-ci.org/cs207G6/cs207-FinalProject.svg?branch=master&maxAge=0)](https://travis-ci.org/cs207G6/cs207-FinalProject.svg?branch=master&maxAge=0)

[![Coverage Status](https://coveralls.io/repos/github/cs207G6/cs207-FinalProject/badge.svg?branch=master)](https://coveralls.io/github/cs207G6/cs207-FinalProject?branch=master&maxAge=0)

# CS207 Group 6
## Overview: 
- The **chemkin** library aims to help solve chemical kinetics problems, such as calculating progress rates and reaction rates.
- We also provide **UI** and **Web API** options to generate plots and return calculation results. To access these features, simply follow the instructions below.

## Documentation Link: 
- https://github.com/cs207G6/cs207-FinalProject/blob/master/documentation/cs207-model-doc.pdf

## Installation Steps:
- Clone the repository: **git clone https://github.com/cs207G6/cs207-FinalProject.git** to your desired directory.
- Change working directory to the root directory of the cloned repository
- Install using **pip install .** or **python setup.py install**
- If desired, run tests using **python setup.py test**
- To start the web UI, type: **python -c "import chemkin.webserver; chemkin.webserver.WebServer(8080).start()"**. Copy the link **http://127.0.0.1:8080/** to your web browser.