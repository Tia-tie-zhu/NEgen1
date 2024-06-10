
# Welcome to  NEgen1  

![Version](https://img.shields.io/badge/version-1.0-blue.svg?cacheSeconds=2592000)
![Prerequisite](https://img.shields.io/badge/c++-8.1.0-blue.svg)

### NEgen1 can describe a range of 2D COFs that follow the nucleation-elongation mechanism

## Prerequisites

cpp >= 8.1.0

## Contents

```sh
NEgen1.cpp # main function
cof_parameter.h # defined the variables in NEgen1.cpp
cof_function.h # defined the expressions of induction time, nucleation rate and growth rate established through symbolic regression
```

## Synthesis condition parameters

Users define parameters in `cof_parameter.h`

```sh
-------------------------------------------
User-defined parameters
-------------------------------------------

int addition_number = 5;# number of addition
int addition_interval = 1200; # time interval s
double C_hhtp_add = 0.005; # monomer concentration mol/L
double V_add[] = {0.02,0.06,0.1777776,0.1066668,0.0355556}; # volume of addition L

```

## Author

👤 **Jiaxin Tian**

👤 **Haoyuan Li**

## Contributing

Contributions, issues and feature requests are welcome!


