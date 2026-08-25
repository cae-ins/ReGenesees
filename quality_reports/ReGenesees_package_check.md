# Rapport de validation — ReGenesees 2.5.1 CAE-INS

Date : 26 août 2026
Amont comparé : `DiegoZardetto/ReGenesees@1432c55be5ed104a44023c35b504443c46d5d584`

## Verdict

La version CAE-INS apporte des améliorations numériques et de sûreté réelles,
sans différence statistique détectée sur les scénarios de compatibilité. Elle
est apte à être publiée sur GitHub comme version locale CAE-INS.

## Résultats numériques

- Calage : 10 scénarios publics comparés à l'amont ; codes retour identiques ;
  écart relatif maximal sur les poids `8,203162e-16`.
- `svyDelta` : 10 scénarios élémentaires, par grappes et par domaines ; écart
  numérique maximal `0` ; structures de sortie identiques.
- Solveur isolé : accélération observée de `2,14x` à `2,75x` selon le passage.
- `svyDelta` multi-domaines : accélération observée de `1,18x` à `1,62x`.
- Raking complet : écart maximal de poids `2,220446e-16`, erreur relative de
  contrainte `4,265511e-12`. Aucun gain de vitesse n'est revendiqué.
- Cas extrême : le Newton amorti converge avec 22 réductions de pas, tandis que
  l'algorithme amont non amorti ne converge pas.

## Contrôles de package

- `devtools::test()` : réussi, 37 attentes, 0 échec.
- `R CMD build` : réussi ; archive `ReGenesees_2.5.1.tar.gz` créée.
- `R CMD check --as-cran` : installation, chargement, déchargement, analyse du
  code, documentation, exemples, tests et manuels réussis.
- Résultat local final : 1 WARNING environnemental (`qpdf` absent pour contrôler
  la compression des PDF) et 4 NOTE non fonctionnelles (nouvelle soumission,
  deux exemples dépassant 5 secondes sur cette machine, outils HTML
  `tidy`/`V8` absents, fichier temporaire MiKTeX). Aucun ERROR, aucune anomalie
  du code R et aucun lien CAE-INS invalide après création du dépôt.

## Revue indépendante et corrections

La revue source indépendante avait notamment signalé des écritures dans
`.GlobalEnv`, une modification globale des contrastes, le déverrouillage du
namespace GVF, un helper assimilable à tort à une méthode S3, un import
`memory.limit` obsolète et un démarrage fondé sur `sink()`.

Ces points ont été corrigés :

- diagnostics attachés aux objets retournés ou aux conditions d'erreur ;
- contrastes ReGenesees appliqués localement aux matrices de modèle ;
- état GVF conservé dans un environnement mutable sans `unlockBinding()` ;
- helper renommé `.svyDelta_by()` ;
- imports `MASS` et `memory.limit` supprimés ;
- message de démarrage produit avec `capture.output()`.

La migration complète de la documentation historique vers roxygen2 n'a pas été
entreprise : elle n'améliore pas les algorithmes, modifierait mécaniquement une
grande surface du package et n'est pas nécessaire au contrôle réussi de cette
version. Les pages Rd concernées par les changements d'API ont été mises à jour
manuellement, conformément à la structure du projet amont.

## Reproduction

```powershell
Rscript -e "devtools::test(reporter='summary', stop_on_failure=TRUE)"
Rscript tools/benchmarks/verify-upstream-compatibility.R compare <baseline.rds>
Rscript tools/benchmarks/verify-delta-compatibility.R compare <baseline.rds>
Rscript tools/benchmarks/benchmark-calibration.R
R CMD build ReGenesees-CAE-INS
R CMD check --as-cran ReGenesees_2.5.1.tar.gz
```
