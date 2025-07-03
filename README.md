# ADMSFI
ADMSFI : Anomaly Detection Based on Multi-Sequence Fuzzy Feature Interaction

## Usage

You can run `Demo_ADMSFI.m`

```matlab
clear all
clc

format short
load Example.mat

Dataori=Example;

trandata=Dataori;
trandata(:,2:3)=normalize(trandata(:,2:3),'range');
sigma=1;

out_scores=ADMSFI(trandata,sigma)
```
You can get outputs as follows:

```matlab
out_scores =

    0.8132
    0.8024
    0.8132
    0.7995
    0.8995
    0.8048
```
