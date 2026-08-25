[![Référencé dans Awesome Official Statistics](https://awesome.re/mentioned-badge.svg)](https://github.com/SNStatComp/awesome-official-statistics-software)

# ReGenesees <img src="ReGenesees_LOGO_small.png" align="right" alt="Logo de ReGenesees" />

**ReGenesees (R Evolved Generalized Software for Sampling Estimates and Errors
in Surveys)** est un package R consacré à l'analyse *design-based* et
*model-assisted* des enquêtes par sondage complexes.

> **Ce dépôt est notre version locale CAE-INS de ReGenesees.** Il s'agit d'un
> fork maintenu par CAE-INS, dérivé de
> [DiegoZardetto/ReGenesees](https://github.com/DiegoZardetto/ReGenesees).
> Le 25 août 2026, CAE-INS a amélioré les solveurs de calage, vectorisé des
> opérations de `svyDelta` et ajouté des tests, benchmarks et contrôles
> automatisés. Ce fork n'est pas une version officielle de l'Istat et n'est pas
> approuvé par les auteurs d'origine. Voir [NOTICE.md](NOTICE.md) pour
> l'attribution détaillée et le périmètre des modifications.

ReGenesees est l'outil de référence de l'**Istat (Institut national de
statistique italien)** pour le calage, l'estimation et l'évaluation des erreurs
d'échantillonnage.

## Installation

Installer directement la dernière version CAE-INS depuis GitHub :

```r
# Méthode recommandée
install.packages("pak")
pak::pak("cae-ins/ReGenesees")

# Autre méthode
install.packages("remotes")
remotes::install_github("cae-ins/ReGenesees")
```

Les versions étiquetées sont également disponibles sous forme de package source
`.tar.gz` et de binaire Windows `.zip` sur la page
[GitHub Releases](https://github.com/cae-ins/ReGenesees/releases).

Le fork CAE-INS conserve le nom R `ReGenesees` afin de préserver la
compatibilité. Il remplace donc toute autre version installée de `ReGenesees`
dans une même bibliothèque R.

## Améliorations vérifiées

La version 2.5.1 apporte notamment :

- un solveur hybride utilisant Cholesky dans les systèmes réguliers et une SVD
  tronquée dans les systèmes singuliers ou mal conditionnés ;
- un Newton amorti pour éviter les divergences et débordements numériques du
  raking et des calages bornés ;
- la suppression de calculs redondants dans les itérations de calage ;
- une exécution plus efficace de `svyDelta` par domaine : indices calculés une
  seule fois, agrégation des contributions PSU avec `rowsum()` et calcul à la
  demande des mesures de variabilité ;
- l'extension du solveur hybride aux chemins historiques de calage et aux
  linéarisations de coefficients de régression ;
- un chargement plus sûr, sans manipulation fragile de `sink()` ;
- la suppression des écritures implicites dans `.GlobalEnv`, du déverrouillage
  du namespace et des modifications globales de contrastes ;
- une suite de tests numériques et des comparaisons systématiques avec l'amont.

Les preuves reproductibles donnent les résultats suivants :

| Composant | Compatibilité avec l'amont | Gain mesuré |
|---|---:|---:|
| Solveur linéaire isolé | écart maximal `2,22e-16` | de **2,14× à 2,75×** selon le passage |
| Raking complet | erreur de contrainte `4,27e-12` | résultats équivalents et convergence renforcée |
| `svyDelta` multi-domaines | écart numérique **nul** sur 10 cas | de **1,18× à 1,62×** selon le passage |

Le Newton amorti converge aussi sur un cas extrême de raking où l'algorithme
amont non amorti échoue. Comme les chronométrages dépendent de la charge de la
machine, les intervalles ci-dessus reprennent plusieurs passages ; aucun gain de
vitesse n'est revendiqué pour le raking complet. Les chiffres exacts, les jeux
de comparaison et les scripts reproductibles se trouvent dans
`quality_reports/` et `tools/benchmarks/`.

## Versions officielles amont

Les versions officielles sont diffusées sur le
[site de l'Istat](https://www.istat.it/en/classifications-and-tools/methods-and-software-of-the-statistical-process/process-phase/weighting-estimation-and-sampling-error-evaluation/regenesees/)
et sur la
[plateforme Joinup de la Commission européenne](https://joinup.ec.europa.eu/solution/regenesees-system/releases),
qui conserve également les versions antérieures.

## Documentation

Le site `pkgdown` du projet amont est disponible à l'adresse suivante :
<https://diegozardetto.github.io/ReGenesees/>.

## Citation

Zardetto, D. (2015). « ReGenesees: An Advanced R System for Calibration,
Estimation and Sampling Error Assessment in Complex Sample Surveys ».
*Journal of Official Statistics*, 31(2), 177–203.
<https://reference-global.com/article/10.1515/jos-2015-0013>.

## Interface graphique

Le package compagnon
[ReGenesees.GUI](https://github.com/DiegoZardetto/ReGenesees.GUI) fournit une
interface graphique Tcl/Tk pour utiliser ReGenesees à la souris.

## Soutiens du projet amont

Le projet ReGenesees a été conçu à l'Istat à la fin de l'année 2006 ; l'Istat
en est depuis le principal soutien. Depuis avril 2021, le développement, la
maintenance et l'assistance du projet amont sont également soutenus activement
par la Banque mondiale.

## Licence et attribution

Le code est distribué sous la **Licence publique de l'Union européenne (EUPL)**.
Le texte intégral figure dans [LICENSE](LICENSE). Les auteurs et copyrights
d'origine sont conservés ; les modifications CAE-INS sont identifiées dans
[NOTICE.md](NOTICE.md), [NEWS.md](NEWS.md) et l'historique Git.
