#!/bin/bash

if [ -x /opt/local/bin/convert ]
then
    CONVERT=/opt/local/bin/convert
elif [ -x /usr/local/bin/convert ]
then
    CONVERT=/usr/local/bin/convert
else
    exit 2
fi

if [ -d ./images ]
then
    IMAGEDIR=./images
else
    exit 3
fi

if [ -x ${CONVERT} ]
then
	${CONVERT} -delay 20 ${IMAGEDIR}/Pressure*png AnimPressure.gif
	${CONVERT} -delay 20 ${IMAGEDIR}/Stream*png   AnimStream.gif
	${CONVERT} -delay 20 ${IMAGEDIR}/Forces*png   AnimForces.gif
	${CONVERT} -delay 20 ${IMAGEDIR}/Velocity*png AnimVelocity.gif
	${CONVERT} -delay 20 ${IMAGEDIR}/ParticleVelocity*png AnimParticles.gif
	${CONVERT} -delay 20 ${IMAGEDIR}/ParticleTrace*png AnimParticleTrace.gif
else
	echo "Install the ImageMagick package and adapt this script to locate convert"
	exit 1
fi

