#!/bin/bash

# Modify these parameters
# -----------------------------------------------------------------------------
trace_particles=true
save_sorted_files=true
load_tracer_meta=false
run_name=mime1836_sigmaic016_bg00
runpath=/lustre/scratch4/turquoise/.mdt0/xiaocan/reconnection/power_law_index/mime1836/$run_name
filepath=$runpath/tracer/tracer1/
particle=electron
tstep_min=0
tstep_max=69683
tstep_interval=17
is_recreate=0 # recreate a file?
nsteps=100
reduced_tracer=0
mpi_size=36
ux_index=6
q_index=13       # particle tag index in the HDF5 file
energy_index=14  # >= number of datasets in the HDF5 file
ntraj=1000       # number of trajectories to search for
ratio_emax=1     # maximum energy / starting energy
echo "Maximum time step:" $tstep_max
echo "Time interval:" $tinterval

tstep=$tstep_max
input_file=tracers.h5p
meta_file=tracers.h5p
energy_sorted_file=${particle}_tracer_energy_sorted.h5p
qtag_sorted_file=${particle}_tracer_qtag_sorted.h5p
group_name=/Step#$tstep
group_name_output=/Step#$tstep
subgroup_name=/Step#$tstep/${particle}_tracer
meta_group_name=/Step#$tstep/grid_metadata
single_h5=1  # tracers for all species + metadata are saved in a single file
single_group=1  # nsteps of tracer data are saved in the same group
# -----------------------------------------------------------------------------

if [ "$trace_particles" = true ] ; then
    additional_flags="-q"
fi

if [ "$save_sorted_files" = false ] ; then
    additional_flags="$additional_flags -w"
fi

if [ "$load_tracer_meta" = true ] ; then
    additional_flags="-r $additional_flags"
fi

echo "Additional flags: " $additional_flags

cd ../

# Sort the tracer particle at $tstep ($tstep_max in default) at the by particle
# energies. This will generate a new file $energy_sorted_file in the directory
# ($fpath below). Particles at the beginning of $energy_sorted_file have the
# Lowest energy. Particles at the end of of $energy_sorted_file have the highest
# energy.
tinterval_file=$(($nsteps * $tstep_interval))
tstep_file=$(( $tstep / $tinterval_file ))
tstep_file=$(( $tstep_file * $tinterval_file))
fpath=$filepath/T.$tstep_file
rm $fpath/$energy_sorted_file
srun -n $mpi_size \
./h5group-sorter -f $input_file \
                 -o $energy_sorted_file \
                 -g $group_name \
                 -m $meta_file \
                 -k $energy_index \
                 -a attribute \
                 -u $ux_index \
                 --tmin=$tstep_min \
                 --tmax=$tstep_max \
                 --tstep=$tstep \
                 --tinterval=$tstep_interval \
                 --filepath=$filepath \
                 --species=$particle \
                 --is_recreate=$is_recreate \
                 --nsteps=$nsteps \
                 --reduced_tracer=$reduced_tracer \
                 --single_h5=$single_h5 \
                 --single_group=$single_group \
                 --subgroup_name=$subgroup_name \
                 --meta_group_name=$meta_group_name \
                 --group_name_output=$group_name_output

# Sort tracer particles at all time steps by particle tags (q dataset in the
# *.h5p files) and save the sorted tracers into $qtag_sorted_file.
# At the same time, we will get some particle trajectories.
# 1. We will select some high-energy particles from $energy_sorted_file based on
#    $ratio_emax. Assuming the highest-energy particle has energy emax, we will
#    search for particles with energies closest to but smaller than emax/ratio_emax.
# 2. We will keep tracking the these selected tracer particles as we sort through
#    the tracer particles at all time steps.
# 3. We will save the tracer trajectory data into $filename_traj below. Each tracer
#    particle in $filename_traj will occupy one group.

# filename for particle trajectory data
filename_traj=data/${particle}s_ntraj${ntraj}_${ratio_emax}emax.h5p

srun -n $mpi_size \
./h5group-sorter -f $input_file \
                 -o $qtag_sorted_file \
                 -g $group_name \
                 -m $meta_file \
                 -k $q_index \
                 -a attribute \
                 -u $ux_index \
                 -p $additional_flags \
                 --tmin=$tstep_min \
                 --tmax=$tstep_max \
                 --tstep=$tstep \
                 --tinterval=$tstep_interval \
                 --filepath=$filepath \
                 --species=$particle \
                 --filename_traj=$filename_traj \
                 --nptl_traj=$ntraj \
                 --ratio_emax=$ratio_emax \
                 --is_recreate=$is_recreate \
                 --nsteps=$nsteps \
                 --reduced_tracer=$reduced_tracer \
                 --single_h5=$single_h5 \
                 --single_group=$single_group \
                 --subgroup_name=$subgroup_name \
                 --meta_group_name=$meta_group_name \
                 --group_name_output=$group_name_output

data_dir=data/power_law_index/$run_name
mkdir -p $data_dir
mv $filename_traj $data_dir

cd config
