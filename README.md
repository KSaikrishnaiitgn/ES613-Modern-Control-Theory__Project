# Modeling and Control of a Soft Robotic Gripper  
### Using LQR with Observer and Disturbance Handling

---

## 📌 Overview
This project focuses on the modeling, analysis, and control of a soft robotic gripper using modern control techniques. A nonlinear system is first linearized, followed by system analysis and controller design using Linear Quadratic Regulator (LQR). An observer is incorporated for state estimation, and disturbance effects are studied to evaluate robustness.

---

## ⚙️ Methodology
The workflow of the project is as follows:

1. **System Linearization**  
   The nonlinear dynamics are linearized using Jacobian-based methods.

2. **System Analysis**  
   The linearized model is analyzed for stability, controllability, and observability.

3. **LQR Controller,Observer Design**  
   An optimal state-feedback controller is designed with output weighting.
   A state observer is developed to estimate unmeasured states.

4. **Disturbance Handling**  
   The system is evaluated under disturbances to test robustness.

5. **Noise Handling**
    The system is evaluated under noise to test robustness.
---

## 📁 File Structure
1. **01_Jacobian_linearization.m**
2. **02_System_Analysis.m** 
3. **03_lQR_Obs_design.m**
4. **04_Gripper_Autodisturbance.m**
5. **05_Gripper_Outputweighted.m**
6. **README.md**
## ▶️ How to Run
1. Open MATLAB and navigate to the project folder  
2. Run the files in the following order:
1. **01_Jacobian_linearization.m**
2. **02_System_Analysis.m** 
3. **03_lQR_Obs_design.m**
4. **04_Gripper_Autodisturbance.m**
5. **05_Gripper_Outputweighted.m**

## 🧠 Key Concepts Used

- Nonlinear System Linearization (Jacobian Method)  
- State-Space Representation  
- Controllability and Observability  
- Linear Quadratic Regulator (LQR)  
- State Observer Design  
- Disturbance Rejection  
- Performance for Noise application
---

## 📊 Expected Outcomes

- Linearized state-space model  
- Stability and controllability insights  
- Optimal feedback gain matrix  
- Accurate state estimation using observer  
- Improved robustness under disturbances  

---

## 🛠️ Requirements

- MATLAB (R2020 or later recommended)  
- Control System Toolbox  

---

## 👨‍💻 Author

K. Saikrishna  

---

## 📚 Notes

This project is developed as part of coursework in **Modern Control Theory**. The implementation emphasizes clarity, modularity, and reproducibility.

---
