# Esto hace un ambiente containerizado
# usaremos mamba para recrearlo en un contenedor
# puede demorar :/ 

FROM python

# Define la carpeta de trabajo
WORKDIR /usr/src/app

# Instala cosas en el base
RUN pip install --no-cache-dir \
    pandas \
    numpy \
    scikit-learn \
    cobra \
    networkx \
    matplotlib \
    ray \
    numba

# Corre los comandos y cosas
# CMD "ray start --head && \
#    python ./src/step02_centralidades_computo.py"

# Copia todo dentro del entorno del contenedor
COPY . .