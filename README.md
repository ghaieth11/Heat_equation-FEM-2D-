# Résolution de l'Équation de la Chaleur en 1D par la Méthode des Éléments Finis  

## Description  
Ce projet vise à résoudre l'équation de la chaleur en une dimension à l'aide de la **méthode des éléments finis (EF)**. L'équation étudiée est la suivante :  


$- u''(x) = f(x), \quad 0 < x < L,$


avec des **conditions aux limites de Dirichlet homogènes** :  

$$
u(0) = 0, \quad u(L) = 0.
$$

L'objectif est d'analyser le problème théoriquement et de l'implémenter numériquement en **Julia**.

## Structure du Projet  

Le travail est divisé en deux parties principales :  

### 1. Analyse Théorique  
- Formulation forte et formulation faible du problème.  
- Définition de l’espace d’approximation de dimension finie.  
- Construction de la **matrice de rigidité** et du **système linéaire** associé.  
- Étude de deux cas de source :  
  - $f(x) = 1 $  
  - $f(x) = \sin(x)$  

### 2. Résolution Numérique & Implémentation  
- Implémentation de la méthode des éléments finis en **Julia**.  
- Comparaison des solutions :  
  - Solution exacte.  
  - Solution obtenue par l'implémentation numérique.  
  - Solution obtenue par calcul manuel.  
- Étude de l’erreur d’approximation, de la **stabilité** et de la **convergence** de la solution numérique.  
- Discussion des améliorations possibles pour optimiser la précision et l’efficacité de la méthode.

## Prérequis  
Pour exécuter ce projet, vous aurez besoin de Julia. Vous pouvez l’installer via :  

```bash
# Sous Linux/macOS
curl -sSL https://julialang.org/downloads/ | bash

# Sous Windows, téléchargez et installez Julia depuis https://julialang.org/downloads/
```
### Installation et Exécution
```bash
git clone https://github.com/votre-utilisateur/votre-repo.git
cd votre-repo
julia main.jl
```

## Auteur

Ce projet a été réalisé par **ALOUI Ghaieth**, étudiant en **Mathématiques Appliquées et Modélisation** à **POLYTECH Nice Sophia** – **Université Côte d'Azur**, dans le cadre d'un TP sur la méthode des éléments finies **semestre S8**.

