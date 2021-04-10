| The dividend and divisor
dividend: dc #999
divisor: dc #77

| Allocate memory for results
quotient: ds 1
remainder: ds 1

| Allocate data registers in the following manner:
|  cD0 - divisor
|  cD1 - quotient
|  cD2 - remainder
move divisor, D0
move #0, D1
move dividend, D2

| Perform the logic loop
div: inc D1
sub D0, D2
bgt div
beq fin
add D0, D2
dec D1

| We store the results
fin: move D1, quotient
move D2, remainder