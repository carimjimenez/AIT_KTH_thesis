# AIT_KTH_thesis
VSC with Grid-Forming strategies (Droop-Control)

VSC_3phVoltageSource   is a basic model, which includes VSC and 3ph Voltage source

File VSC_GridFroming_Thesis_Script_withSynchMach.m constains all the data of the electric network, parameters of the Inverter and its controllers. This file can be applied to any SIMULINK model

VSC_3phVoltageSource_VirtImp.   a VSC with an implemented Virtual impedance block, which helps to reduce the current magnitude in the inverter during a fault condition (3ph Short-circuit and 1ph-Ground) and a 3ph Voltage source

VSC_VirtualImp_SyncMach  a VSC with an implemented Virtual impedance block  and a Synchronous machine (SM). Uncomment Lines (234-247) in the script to run this file

VSC_VirtualImp_CurrentLimit_SyncMach    a VSC with an implemented Virtual impedance block and current reference limitation  and a Synchronous machine. Uncomment Lines (251-265) in the script to run this file

IEEE9_BUS_withVSC  an IEEE 9-bus system is modelled with a SM and two VSC-GF but one has an adittional Virtual Impedance block
Uncomment Lines (267-285) in the script to run this file
