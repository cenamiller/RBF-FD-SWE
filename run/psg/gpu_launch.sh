#!/bin/bash

export ACC_DEVICE_NUM=$(($OMPI_COMM_WORLD_LOCAL_RANK%$num_tasks))
$*
