# Attribution and modification notice

ReGenesees was created and is maintained upstream by Diego Zardetto. It is the
standard tool for calibration, estimation and sampling-error assessment at the
Italian National Institute of Statistics (Istat). Parts of its survey-analysis
code incorporate modified code derived from early versions of Thomas Lumley's
`survey` package.

Upstream source: <https://github.com/DiegoZardetto/ReGenesees>

CAE-INS created this derivative from upstream commit
`1432c55be5ed104a44023c35b504443c46d5d584` and modified it on 25 August 2026.
The changes improve the numerical solver used for calibration, add damped
Newton iteration, tests, benchmarks, continuous integration and release
automation. Detailed changes are recorded in `NEWS.md` and in the Git history.

This derivative remains licensed under the European Union Public Licence
(EUPL), whose full text is in `LICENSE`. It is not an official Istat release and
is not endorsed by Diego Zardetto, Thomas Lumley, Istat or the World Bank.
