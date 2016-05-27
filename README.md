# isingmodel

To output all matrix:
  - Compile file:
  ```
  $ make picture=1
  ```
  - Output data:
  ```
  $ ./2dising.sh
  ```
  (If it doesn't work, type in:
  ```
  $ chmod +x ./2dising.sh
  ```
  ).

To output energy, heat capacity, magnetization, magnetization susceptibility:
  - Compile file:
  ```
  $ make
  ```
  - Output data:
  ```
  $ ./2dising.sh
  ```

All data file will be in ./data

To clean all data:
```sh
$ make cleandata
```
How to access Titan box:
```sh
$ ssh username@ieng6.ucsd.edu
$ ssh username@igpu6-210.ucsd.edu
```

How to run the cuda code:
  - Compile file:
  ```
  $ make picture=1
  ```
  - Output data:
  ```
  $ ./2dising.sh
  ```
  (If it doesn't work, type in:
  ```
  $ chmod +x ./2dising.sh
  ```
  ).
