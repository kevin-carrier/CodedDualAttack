################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/CodedDualAttack.c 

C_DEPS += \
./src/CodedDualAttack.d 

OBJS += \
./src/CodedDualAttack.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I/opt/homebrew/Cellar/mpfr/4.1.0/include/ -I/opt/homebrew/Cellar/flint/2.9.0/include/flint/ -I/opt/homebrew/Cellar/gmp/6.2.1_1/include/ -I/opt/homebrew/Cellar/flint/2.9.0/include/flint/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src

clean-src:
	-$(RM) ./src/CodedDualAttack.d ./src/CodedDualAttack.o

.PHONY: clean-src

