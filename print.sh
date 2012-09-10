#!/bin/sh

ROOT="/home/jims/ascidian/solver"

files=""
files="$files RunSims.m"
files="$files 0params/GetParams.m"
files="$files 0params/GetBotry.m"
files="$files 0params/GetKinematics.m"
files="$files solver.m"
files="$files 0attrib/bodyMass.m"
files="$files 0attrib/bodyVolume.m"
files="$files 0attrib/trunkVolume.m"
files="$files 0attrib/inertiaTensor.m"
files="$files 0trans/BodyToInertial.m"
files="$files 0trans/InertialToBody.m"
files="$files 0forces/CalcTailForce.m"
files="$files 0forces/CalcTrunkForce.m"
files="$files 0forces/CalcBouyForce.m"
files="$files 0forces/CalcGravForce.m"

a2ps_opts=""
a2ps_opts="$a2ps_opts --line-numbers=10"
a2ps_opts="$a2ps_opts --font-size=8"
a2ps_opts="$a2ps_opts --columns=1"
a2ps_opts="$a2ps_opts --rows=1"
a2ps_opts="$a2ps_opts -R"
a2ps_opts="$a2ps_opts -b"
a2ps_opts="$a2ps_opts -o -"

for i in $files
do
   a2ps $a2ps_opts $i | lpr
done
