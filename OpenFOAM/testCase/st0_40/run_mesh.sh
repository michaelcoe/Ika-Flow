#!/bin/bash -e

cp -rf ../mesh/backGround/constant/polyMesh ./constant/

rm -rf 0
cp -rf 0_org 0

setFields

setFields -dict system/setFieldsDict_init