"use strict";

function _typeof(obj) { "@babel/helpers - typeof"; if (typeof Symbol === "function" && typeof Symbol.iterator === "symbol") { _typeof = function _typeof(obj) { return typeof obj; }; } else { _typeof = function _typeof(obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; }; } return _typeof(obj); }

require("core-js/modules/es.array.slice");

require("core-js/modules/es.array.some");

exports.__esModule = true;
exports.operationResult = operationResult;
exports.SHORT_NAME_FOR_COMPONENT = exports.OPERATION_ID = void 0;

var C = _interopRequireWildcard(require("../../../i18n/constants"));

var _logicalOperationRegisterer = require("../logicalOperationRegisterer");

function _getRequireWildcardCache() { if (typeof WeakMap !== "function") return null; var cache = new WeakMap(); _getRequireWildcardCache = function _getRequireWildcardCache() { return cache; }; return cache; }

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } if (obj === null || _typeof(obj) !== "object" && typeof obj !== "function") { return { default: obj }; } var cache = _getRequireWildcardCache(); if (cache && cache.has(obj)) { return cache.get(obj); } var newObj = {}; var hasPropertyDescriptor = Object.defineProperty && Object.getOwnPropertyDescriptor; for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) { var desc = hasPropertyDescriptor ? Object.getOwnPropertyDescriptor(obj, key) : null; if (desc && (desc.get || desc.set)) { Object.defineProperty(newObj, key, desc); } else { newObj[key] = obj[key]; } } } newObj.default = obj; if (cache) { cache.set(obj, newObj); } return newObj; }

var OPERATION_ID = 'disjunctionWithExtraCondition';
exports.OPERATION_ID = OPERATION_ID;
var SHORT_NAME_FOR_COMPONENT = C.FILTERS_LABELS_DISJUNCTION; // ((p OR q OR w OR x OR...) AND z) === TRUE?

exports.SHORT_NAME_FOR_COMPONENT = SHORT_NAME_FOR_COMPONENT;

function operationResult(conditions, value) {
  if (conditions.length < 3) {
    throw Error('Operation doesn\'t work on less then three conditions.');
  }

  return conditions.slice(0, conditions.length - 1).some(function (condition) {
    return condition.func(value);
  }) && conditions[conditions.length - 1].func(value);
}

(0, _logicalOperationRegisterer.registerOperation)(OPERATION_ID, SHORT_NAME_FOR_COMPONENT, operationResult);