# Assembleur basé sur les graphes de Debruijn

## Introduction

Ce projet implémente un assembleur génomique simplifié basé sur un graphe de De Bruijn, une méthode largement utilisée dans les pipelines d’assemblage modernes (ex. Velvet, SPAdes).

L’objectif est de reconstruire le génome d’un virus (Enterovirus A71) à partir de reads courts (Illumina).
Le pipeline permet de :
- générer les k-mers et construire le graphe de De Bruijn
- nettoyer le graphe (tips, bulles)
- extraire des contigs
- comparer ces contigs au génome de référence via BLAST

Ce TP illustre le fonctionnement interne des assembleurs modernes, la gestion des erreurs de séquençage et la reconstruction de séquences à partir de fragments.

## Prerequis
Vous devez avoir Python 3.10+ installé,
```
python3 --version
```

Ainsi que le gestionnaire d'environnement 'uv' installé.
```
curl -LsSf https://astral.sh/uv/install.sh | sh
uv --version
```
Ainsi que le logiciel BLAST+ soit installé
```
# Linux
sudo apt install ncbi-blast+
# macOS (Homebrew)
brew install blast              
```

# 


Veuillez cloner le projet et aller au repertoire crée: 
```
git clone https://github.com/anaisdlss/tp_debruijn.git
cd tp_debruijn
```

Sychroniser l'environnement
```
uv sync
```
Créez un fichier *résultats* :
```
mkdir resultats
```

Puis exectutez le script :
```
uv run python debruijn/debruijn.py -i data/eva71_plus_perfect.fq -k 22 -o resultats/contigs.fasta -f resultats/graph.png
```
Il est possible de modifier la taille des kmer avec -k


Après cela, créez la base de données BLAST pour la référence ```data/eva71.fna```:
```
makeblastdb -in data/eva71.fna -dbtype nucl
```
Puis comparez vos contigs obtenues dans ```resultats/contigs.fasta``` avec la référence :
```
blastn -query resultats/contigs.fasta -db data/eva71.fna -out resultats/blast_result.txt -outfmt 6
```

Les resultats obtenus se trouvent dans le fichier ```resultats/blast_result.txt```.
