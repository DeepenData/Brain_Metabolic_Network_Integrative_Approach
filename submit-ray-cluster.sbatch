#!/bin/bash
### Este script funciona para cualquier numero de nodos; Ray se encarga de manejarlos
### puedes definir este parametro en la ejecución de la tarea como 'sbatch submit-ray... --nodes=5'
#SBATCH --nodes=4

### Asigna todos los recursos de un nodo al proceso Ray, el cual los maneja internamente
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-task=0
#SBATCH --cpus-per-task=20

# Carga modulos aqui
# TODO: revisar como funciona esto en NLHPC
# module load miniconda/ 

#conda initrm 
conda activate alejandro

echo "$(date) --- Empieza el setup"

redis_password=$(uuidgen)
export redis_password

nodes=$(scontrol show hostnames $SLURM_JOB_NODELIST) # Getting the node names
nodes_array=( $nodes )

node_1=${nodes_array[0]} 
ip=$(srun --nodes=1 --ntasks=1 -w $node_1 hostname) # making redis-address
port=6379
ip_head=$ip:$port
export ip_head
echo "IP Head: $ip_head"

echo "$(date) --- Empezando HEAD en $node_1"
srun --nodes=1 --ntasks=1 -w $node_1 $(pwd)/src/sbatch-ray/start-head.sh $ip $redis_password &
sleep 30

worker_num=$(($SLURM_JOB_NUM_NODES - 1)) #number of nodes other than the head node
for ((  i=1; i<=$worker_num; i++ ))
do
  node_i=${nodes_array[$i]}
  echo "$(date) --- Empezando WORKER $i en $node_i"
  srun --nodes=1 --ntasks=1 -w $node_i $(pwd)/src/sbatch-ray/start-worker.sh $ip_head $redis_password &
  sleep 5
done

#### SECCIÓN QUE SI EJECUTA CODIGO PYTHON Y DEMAS

echo "$(date) --- Status del cluster Ray"
ray status --address $ip_head --redis_password $redis_password

echo "$(date) --- Empieza el proceso python paralelo"

time python -u $(pwd)/src/step02_centralidades_computo.py

echo "$(date) --- Termina el proceso paralelo"

ray stop --verbose

exit
