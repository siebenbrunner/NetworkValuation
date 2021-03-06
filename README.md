# Network Valuation
Routines for the valuation of cross-holdings of financial contracts.

## Background
This repository provides a set of routines that can be used to determine payoffs and valuations for financial contracts held between a set of institutions in a system, in the spirit of models such as [Eisenberg and Noe (2001)] and [Barucca et al. (2016)]. In its current state, it acts as a companion to the paper by [Puig and Siebenbrunner (2018)], but contributions by other researchers are welcomed.

## Usage
### Installation
If you are familiar with git, fork or clone the repository, otherwise use the download link to obtain the source code for the algorithms. Select an appropriate algorithm (see below) and add the corresponding folder to the Matlab path (*note*: adding all folders to the Matlab path will result in name conflicts).

### Algorithms
The algorithms are organized into folders, which act as separately functionable code bases. Selecting the appropriate algorithm for your problem requires some familiarity with the literature. Currently the following algorithms are implemented (in chronological order, which roughly corresponds to an increasing order of complexity):

* [Eisenberg and Noe (2001)]: computes payoffs for one layer of debt contracts without equity participations.
* [Elsinger (2009)]: extends the Eisenberg/Noe-model to include equity participations and several seniority layers of debt.
* [Rogers and Veraart (2013)]: extends the Eisenberg/Noe-model to include liquidation haircuts.
* [Puig and Siebenbrunner (2018)]: extends the Elsinger-model to include bail-ins and contingent convertible debt securities.

### Examples
Each of the folders contains an example file demonstrating the use of the algorithms. As a simple demonstration, using the Eisenberg/Noe algorithm can be done via a simple function call:
```matlab
vecE = [10;10;10];
matL = [0 20 10;5 0 24;10 0 0];
vecPayments = calcPayments(vecE,matL)
vecPayments =
   19.5455
   26.3636
   10.0000
```

## Contributing
We invite other researchers to contribute their algorithms to the repository, so that it can serve as a common resource for all type of network-valuation-related methodologies. If you are familiar with git you may simply fork the repository and create a pull request, otherwise you can get in touch with us under christoph.siebenbrunner@maths.ox.ac.uk about how to share your code. Contributions may use and extend existing code from this repository. They do not necessarily have to be written in Matlab, as each folder acts as a separate code base. Authorship of contributed code for academic/citation purposes will be fully acknowledged.

## License and scholarly attribution
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details. If used in acadmic works, please cite the companion paper [Puig and Siebenbrunner (2018)], in addition to the original paper, depending on which algorithm is used (*note*: we will update the attribution framework if routines are added by other researchers).

## References
Barucca, P., Bardoscia, M., Caccioli, F., D'Errico, M., Visentin, G., Battiston, S., & Caldarelli, G. (2016). Network valuation in financial systems. SSRN Working Paper  
Eisenberg, L., & Noe, T. H. (2001). Systemic risk in financial systems. Management Science, 47(2), 236-249.  
Elsinger, H. (2009). Financial networks, cross holdings, and limited liability. Oesterreichische Nationalbank Working Paper 156  
Rogers, L. C., & Veraart, L. A. (2013). Failure and rescue in an interbank network. Management Science, 59(4), 882-898.  
Puig, M., & Siebenbrunner, C. (2018) Multilayer network valuation with equity conversions. Working Paper

[Eisenberg and Noe (2001)]: https://doi.org/10.1287/mnsc.47.2.236.9835 "Systemic risk in financial systems"
[Elsinger (2009)]: https://ideas.repec.org/p/onb/oenbwp/156.html "Financial networks, cross holdings, and limited liability"
[Rogers and Veraart (2013)]: https://doi.org/10.1287/mnsc.1120.1569 "Failure and Rescue in an Interbank Network"
[Barucca et al. (2016)]: http://dx.doi.org/10.2139/ssrn.2795583 "Network valuation in financial systems"
[Puig and Siebenbrunner (2018)]: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3225943 "Multilayer network valuation with equity conversions"
