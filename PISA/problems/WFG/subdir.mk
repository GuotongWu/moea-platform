################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../variator.cpp \
../variator_internal.cpp \
../variator_user.cpp \
../wfg_interface.cpp 

OBJS += \
./variator.o \
./variator_internal.o \
./variator_user.o \
./wfg_interface.o 

CPP_DEPS += \
./variator.d \
./variator_internal.d \
./variator_user.d \
./wfg_interface.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


