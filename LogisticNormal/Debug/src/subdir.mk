################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/LogisticNormal.cpp 

OBJS += \
./src/LogisticNormal.o 

CPP_DEPS += \
./src/LogisticNormal.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I"C:\Craig\Models\CASAL2\ThirdParty\boost\boost_1_58_0" -I"C:\Craig\projects\2016\DEE2015-02(TrawlSurveySimulation)\Logistic_normal\Libraries\Eigen" -O0 -g3 -c -fmessage-length=0  -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


