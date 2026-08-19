#include <Adafruit_MPU6050.h>
#include <Adafruit_Sensor.h>
#include <AccelStepper.h>
#include <Wire.h>

Adafruit_MPU6050 mpu;

#define dirPin 7
#define stepPin 6
#define motorInterfaceType 1

AccelStepper stepper = AccelStepper(motorInterfaceType, stepPin, dirPin);

float targetPos = -100;
float maxVelocity = 0;
float prevMaxVelocity = 0;
float maxAcceleration = 0;
float prevGyro[4];
float prevTime[4];
float prevAccel[4];
int numBuffer = 0;
float pauseIter = 0;
bool pauseMotor = false;

void setup(void) {
  Serial.begin(115200);

  // Try to initialize!
  if (!mpu.begin()) {
    Serial.println("Failed to find MPU6050 chip");
    while (1) {
      delay(10);
    }
  }

  // wait for ready
  Serial.println(F("\nSend any character to begin DMP programming and demo: "));
  while (Serial.available() && Serial.read());  // empty buffer
  while (!Serial.available());  // wait for data
  while (Serial.available() && Serial.read());  // empty buffer again

  // set accelerometer range to +-8G
  mpu.setAccelerometerRange(MPU6050_RANGE_16_G);

  // set gyro range to +- 500 deg/s
  mpu.setGyroRange(MPU6050_RANGE_1000_DEG);

  // set filter bandwidth to 21 Hz
  mpu.setFilterBandwidth(MPU6050_BAND_21_HZ);

  stepper.setMaxSpeed(40000);
  stepper.setAcceleration(10000);

  stepper.setSpeed(10000);
  delay(100);
}

void loop() {
  /* Get new sensor events with the readings */
  sensors_event_t a, g, temp;
  mpu.getEvent(&a, &g, &temp);

  Serial.print("Sens");
  Serial.print(",");

  /* Print out the values */

  float velN = 0.5;
  float accelN = 0.1;
  float velThresh = 1;      // squat
  //float accelThresh = (0.638 + 0.1)*pow(prevMaxVelocity/10,2)+0.5;  // stand
  float accelThresh = 1;
  float position = 100;

  float accelS = a.acceleration.z+0.3;
  float gyro = g.gyro.x + 0.1;
  float previousGyro = prevGyro[0];
  float time = g.timestamp;
  float timestep = time - prevTime[3];
  float changeGyro = 90*(gyro - prevGyro[1])/(timestep);
  float logicOne = changeGyro/gyro;
  float accel = changeGyro;

  //Serial.print("Var");
  //Serial.print(",");

  //numBuffer = numBuffer-1;

  pauseMotor = false;
  if (timestep > 50){
    pauseMotor = true;
    pauseIter = 1;
  }

  if (pauseIter == 0) {
    pauseMotor = false;
  } else if (pauseIter <= 10) {
    pauseMotor = true;
    pauseIter = pauseIter+1;
  } else {
    pauseIter = 0;
    pauseMotor = false;
  }

  //pauseMotor = true;
    
  //Serial.print("IterPause");
  //Serial.print(",");

  //else if (numBuffer > 0)
    //pauseMotor = true;
  
  //Serial.print(a.acceleration.x);
  //Serial.print(",");
  //Serial.print(timestep/200);
  //Serial.print(",");
  //Serial.print(-2);
  //Serial.print(",");
  //Serial.print(2);
  //Serial.print(",");
  //Serial.print(stepper.currentPosition()/position+0.5);
  //Serial.print(",");
  //Serial.print(accel);
  //Serial.print(",");
  //Serial.print(a.acceleration.z);
  //Serial.print(", ");
  //Serial.print(pauseIter);
  //Serial.print(",");
  //Serial.print(prevMaxVelocity);
  //Serial.print(",");
  //Serial.print(gyro);
  //Serial.print(",");
  //Serial.print(changeGyro);
  //Serial.print(",");
  //Serial.print(accelS);
  //Serial.println(",");
  //Serial.print(numBuffer);
  //Serial.println("");
  //Serial.print(g.gyro.y);
  //Serial.print(",");
  //Serial.print(g.gyro.z);
  //Serial.println("");

  if (pauseMotor){
  } else if (((gyro < -velN) && (prevAccel[2] <= -accelThresh && accel >= -accelThresh)) || ((gyro > velN) && (prevAccel[2] >= accelThresh && accel <= accelThresh))) {
    stepper.runToNewPosition(0);  // stand
    numBuffer = 0;
  } else if (((accel > accelN) && (prevGyro[2] <= -velThresh && gyro >= -velThresh)) || ((accel < -accelN) && (prevGyro[2] >= velThresh && gyro <= velThresh))) {
    stepper.runToNewPosition(-position);  // squat
    numBuffer = 0;
  } else {
  }
  stepper.run();

  //Serial.print("Motor");
  //Serial.print(",");

  if (gyro*prevGyro[0] >= 0) { // If swing is going in one direction
    if (abs(gyro) >= maxVelocity){  // Update max velocity if necessary
        maxVelocity = abs(gyro);
    } else {}
  } else { // If swing changed direction
    prevMaxVelocity = maxVelocity;  // Store previous max velocity
    maxVelocity = 0;                // Reset max velocity
  }

  //Serial.print("MaxVel");
  //Serial.print(",");

  prevGyro[3] = prevGyro[2];
  prevGyro[2] = prevGyro[1];
  prevGyro[1] = prevGyro[0];
  prevGyro[0] = gyro;

  prevAccel[3] = prevAccel[2];
  prevAccel[2] = prevAccel[1];
  prevAccel[1] = prevAccel[0];
  prevAccel[0] = accel;

  prevTime[3] = prevTime[2];
  prevTime[2] = prevTime[1];
  prevTime[1] = prevTime[0];
  prevTime[0] = time;

  Serial.print(prevMaxVelocity);
  Serial.print(",");
  Serial.print(accelThresh);
  Serial.print(",");
  Serial.print(-accelThresh);
  Serial.print(",");
  Serial.print(gyro);
  Serial.print(",");
  Serial.print(accel);
  Serial.println("");
}