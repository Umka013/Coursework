#pragma once

#include <Eigen/Dense>
#include <cmath>

namespace fem {
template <typename T>
class FiniteElement {
 public:
  using Value = T;
  using Size = unsigned long long int;

  Size const static nNodes = 4;
  Size const static nNodeDofs = 2;
  Size const static nElemDofs = nNodes * nNodeDofs;
  Size const static nVoigt = 3;
  Size const static nSub = 3;

  using Vector = Eigen::VectorX<Value>;
  ;

  using Coordinates = Eigen::Vector<Value, nNodeDofs>;
  using Nodes = Eigen::Vector<Coordinates, nNodes>;
  using StiffnessMatrix = Eigen::Matrix<Value, nElemDofs, nElemDofs>;
  using ElasticityMatrix = Eigen::Matrix<Value, nVoigt, nVoigt>;
  using ShapeFunctions = Eigen::Vector<Value, nNodes>;
  using DifferentiationMatrix = Eigen::Matrix<Value, nVoigt, nElemDofs>;
  using DeformationVectors = Eigen::Matrix<Value, nVoigt, nNodes>; 

  void integrateSubfields(
      StiffnessMatrix &sm, Nodes const &nodes, Value const &elasticityModulus,
      Value const &poissonRatio, DifferentiationMatrix &dm) {
    sm.setZero();
    Coordinates lsX, lsY;
    this->findLimits(lsX, lsY, nodes);

    Coordinates od;
    this->calculateOverallDimensions(od, lsX, lsY);
    Value const &a = od(0), &b = od(1);

    ElasticityMatrix em;
    this->calculateElasticityMatrix(em, elasticityModulus, poissonRatio);
    Value const aH = a * 0.5;
    Value const bH = b * 0.5;

    Value xi, eta, weight;

    Value const subA = a / nSub, subB = b / nSub;
    Value const halfSubA = subA * 0.5, halfSubB = subB * 0.5;
    Value subFieldCenterX = -aH + halfSubA, subFieldCenterY = -bH + halfSubB;
    Value xi_glob, eta_glob;

    const Value gaussPoints[2] = {
        (-1.0 / std::sqrt(3.0)) * halfSubA, (1.0 / std::sqrt(3.0)) * halfSubB};
    const Value gaussWeights[2] = {1.0 * halfSubA, 1.0 * halfSubB};

    for (Size subI = 0; subI < nSub; subI++) {    // x
      for (Size subJ = 0; subJ < nSub; subJ++) {  // y
        for (int i = 0; i < nNodeDofs; ++i) {
          for (int j = 0; j < nNodeDofs; ++j) {
            xi = gaussPoints[i];
            eta = gaussPoints[j];
            weight = gaussWeights[i] * gaussWeights[j];

            xi_glob = subFieldCenterX + xi;
            eta_glob = subFieldCenterY + eta;

            this->calculateDifferentiationMatrix(dm, xi_glob, eta_glob, a, b);

            if (subI == 1 and subJ == 1) {
              sm += dm.transpose() * em * 1e-6 * dm * weight;
            } else {
              sm += dm.transpose() * em * dm * weight;
            }
          }
        }
        subFieldCenterY += subB;
      }
      subFieldCenterX += subA;
      subFieldCenterY = -bH + halfSubB;
    }
  }
  void calculateDifferentiationMatrix(
      DifferentiationMatrix &dm, Value const &xi, Value const &eta,
      Value const &a, Value const &b) {
    ShapeFunctions sf_dx, sf_dy;
    this->calculateShapeFunctions_dx(sf_dx, xi, eta, a, b);
    this->calculateShapeFunctions_dy(sf_dy, xi, eta, a, b);

    dm.setZero();
    for (int i = 0; i < nNodes; ++i) {
      dm(0, 2 * i) = sf_dx(i);      // du/dx
      dm(1, 2 * i + 1) = sf_dy(i);  // dv/dy
      dm(2, 2 * i) = sf_dy(i);      // du/dy
      dm(2, 2 * i + 1) = sf_dx(i);  // dv/dx
    }
  }

  DeformationVectors calculateNDS(
      Vector const &displacementVector,
      Value const &elastMod,
      Value const &poissonRat,
      Nodes const &nodes) {
    Coordinates lsX, lsY;
    this->findLimits(lsX, lsY, nodes);
    
    Coordinates od;
    this->calculateOverallDimensions(od, lsX, lsY);
    Value const &a = od(0), &b = od(1);
    
    Value xi_loc = a/2;
    Value eta_loc = b/2;
    
    ElasticityMatrix em;
    calculateElasticityMatrix(em, elastMod, poissonRat);
        
    DifferentiationMatrix dm1;
    DifferentiationMatrix dm2;
    DifferentiationMatrix dm3;
    DifferentiationMatrix dm4;
    calculateDifferentiationMatrix(dm1, -xi_loc, -eta_loc, a, b);
    calculateDifferentiationMatrix(dm2, xi_loc, -eta_loc, a, b);
    calculateDifferentiationMatrix(dm3, xi_loc, eta_loc, a, b);
    calculateDifferentiationMatrix(dm4, -xi_loc, eta_loc, a, b);
    DeformationVectors deformationVectors;

    deformationVectors.setZero();
    deformationVectors.col(0) = dm1 * displacementVector;
    deformationVectors.col(1) = dm2 * displacementVector;
    deformationVectors.col(2) = dm3 * displacementVector;
    deformationVectors.col(3) = dm4 * displacementVector;
    return em * deformationVectors;
  }

    DeformationVectors calculateDifVect(
      Vector const &displacementVector,
      Value const &elastMod,
      Value const &poissonRat,
      Nodes const &nodes) {
    Coordinates lsX, lsY;
    this->findLimits(lsX, lsY, nodes);
    
    Coordinates od;
    this->calculateOverallDimensions(od, lsX, lsY);
    Value const &a = od(0), &b = od(1);
    
    Value xi_loc = a/2;
    Value eta_loc = b/2;
    
    ElasticityMatrix em;
    calculateElasticityMatrix(em, elastMod, poissonRat);
        
    DifferentiationMatrix dm1;
    DifferentiationMatrix dm2;
    DifferentiationMatrix dm3;
    DifferentiationMatrix dm4;
    calculateDifferentiationMatrix(dm1, -xi_loc, -eta_loc, a, b);
    calculateDifferentiationMatrix(dm2, xi_loc, -eta_loc, a, b);
    calculateDifferentiationMatrix(dm3, xi_loc, eta_loc, a, b);
    calculateDifferentiationMatrix(dm4, -xi_loc, eta_loc, a, b);
    DeformationVectors deformationVectors;

    deformationVectors.setZero();
    deformationVectors.col(0) = dm1 * displacementVector;
    deformationVectors.col(1) = dm2 * displacementVector;
    deformationVectors.col(2) = dm3 * displacementVector;
    deformationVectors.col(3) = dm4 * displacementVector;
    return deformationVectors;
  }

 private:
  void calculateShapeFunctions(
      ShapeFunctions &sf, Value const &xi, Value const &eta, Value const &a,
      Value const &b) {
    sf(0) = (a / 2 - xi) * (b / 2 - eta) / (a * b);
    sf(1) = (a / 2 + xi) * (b / 2 - eta) / (a * b);
    sf(2) = (a / 2 + xi) * (b / 2 + eta) / (a * b);
    sf(3) = (a / 2 - xi) * (b / 2 + eta) / (a * b);
  }

  void calculateShapeFunctions_dx(
      ShapeFunctions &sf_dx, Value const &xi, Value const &eta, Value const &a,
      Value const &b) {
    sf_dx(0) = -(b / 2 - eta) / (a * b);
    sf_dx(1) = (b / 2 - eta) / (a * b);
    sf_dx(2) = (b / 2 + eta) / (a * b);
    sf_dx(3) = -(b / 2 + eta) / (a * b);
  }

  void calculateShapeFunctions_dy(
      ShapeFunctions &sf_dy, Value const &xi, Value const &eta, Value const &a,
      Value const &b) {
    sf_dy(0) = -(a / 2 - xi) / (a * b);
    sf_dy(1) = -(a / 2 + xi) / (a * b);
    sf_dy(2) = (a / 2 + xi) / (a * b);
    sf_dy(3) = (a / 2 - xi) / (a * b);
  }

  void calculateElasticityMatrix(
      ElasticityMatrix &em, Value const &elasticityModulus,
      Value const &poissonRatio) {
    Value factor = elasticityModulus / (1 - poissonRatio * poissonRatio);
    em << 1, poissonRatio, 0, poissonRatio, 1, 0, 0, 0, (1 - poissonRatio) / 2;
    em *= factor;
  }

  void findLimits(Coordinates &lsX, Coordinates &lsY, Nodes const &nodes) {
    Value minx = nodes(0)(0), maxx = nodes(0)(0);
    Value miny = nodes(0)(1), maxy = nodes(0)(1);
    for (int i = 0; i < nNodes; ++i) {
      if (nodes(i)(0) < minx) minx = nodes(i)(0);
      if (nodes(i)(0) > maxx) maxx = nodes(i)(0);
      if (nodes(i)(1) < miny) miny = nodes(i)(1);
      if (nodes(i)(1) > maxy) maxy = nodes(i)(1);
    }
    lsX(0) = minx;
    lsX(1) = maxx;
    lsY(0) = miny;
    lsY(1) = maxy;
  }

  void calculateOverallDimensions(
      Coordinates &od, Coordinates const &lsX, Coordinates const &lsY) {
    od(0) = lsX(1) - lsX(0);
    od(1) = lsY(1) - lsY(0);
  }

  void calculateLocalNodes(
      Nodes &locNodes, Coordinates const &lsX, Coordinates const &lsY,
      Nodes const &nodes) {
    Value avgX = (lsX(0) + lsX(1)) / 2;
    Value avgY = (lsY(0) + lsY(1)) / 2;
    for (int i = 0; i < nNodes; ++i) {
      locNodes(i)(0) = nodes(i)(0) - avgX;
      locNodes(i)(1) = nodes(i)(1) - avgY;
    }
  }
};
}  // namespace fem