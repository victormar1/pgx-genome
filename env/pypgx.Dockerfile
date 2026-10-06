# Image du typeur complementaire (etage 4b).
#
# PyPGx etait la seule dependance du module a vivre sur l'hote : une
# installation Python a maintenir sur chaque noeud, que les plateformes
# refusent a juste titre. Les trois autres outils sont des images epinglees par
# empreinte ; celui-ci le devient.
#
# Deux contraintes dictent le contenu :
#   - python 3.10, et non 3.12 : sur 3.12 le typeur echoue sur un appel pandas
#     retire de la bibliotheque ;
#   - un environnement d'execution Java, parce que le phasage statistique passe
#     par beagle, dont l'archive est embarquee dans le paquet PyPGx.
#
# La ressource de donnees (PYPGX_BUNDLE) n'est pas copiee ici : elle pese
# plusieurs gigaoctets, elle evolue independamment du code, et elle se monte.
#
# Construction :
#   docker build -f env/pypgx.Dockerfile -t pgx-genome/pypgx:0.27.0 env
# Puis relever l'empreinte immuable et la consigner dans env/conteneurs.txt.

FROM python:3.10-slim-bookworm

RUN apt-get update \
 && apt-get install -y --no-install-recommends default-jre-headless \
 && rm -rf /var/lib/apt/lists/*

# Version figee : une montee de version du typeur change les resultats et doit
# etre remesuree (voir doc/VALIDATION.md, section 5).
#
# Une dependance transitive n'est distribuee qu'en source et reclame un
# compilateur. Il est installe et retire dans la meme couche, sans quoi
# l'image porterait une chaine de compilation qu'elle n'utilise jamais.
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential \
 && pip install --no-cache-dir pypgx==0.27.0 \
 && apt-get purge -y --auto-remove build-essential \
 && rm -rf /var/lib/apt/lists/*

# Verification a la construction : le typeur repond, et son archive de phasage
# est bien presente. Un echec ici vaut mieux qu'un echec au premier genome.
RUN pypgx --version \
 && python -c "import glob, os, pypgx; \
d = os.path.dirname(pypgx.__file__); \
j = glob.glob(os.path.join(d, 'api', 'beagle*.jar')); \
assert j, 'archive de phasage absente'; \
print('phasage :', os.path.basename(j[0]))"

ENTRYPOINT ["pypgx"]
