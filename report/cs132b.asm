| The number to check
num: dc #23

| Allocate memory for the result
result: ds 1

| Allocate data registers in the following manner:
|  cD0 - dividend/remainder
|  cD1 - divisor
move num, D0
move #2, D1

| Perform the logic loop
div: sub D1, D0
bgt div
beq fin
inc D1
move num, D0
sub D1, D0
beq prime
move num, D0
jmp div

| Set Prime
prime: move #1, D1

| We store the result
fin: move D1, result
