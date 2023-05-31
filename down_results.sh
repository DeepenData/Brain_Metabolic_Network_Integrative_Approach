#!/bin/bash

ray rsync-down aws-ray-cluster.yml /home/ubuntu/centralities.log aws-downloads/centralities.log
ray rsync-down aws-ray-cluster.yml /home/ubuntu/centralities.pickle aws-downloads/centralities.pickle