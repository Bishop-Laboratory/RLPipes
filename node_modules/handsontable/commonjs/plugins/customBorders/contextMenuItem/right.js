"use strict";

function _typeof(obj) { "@babel/helpers - typeof"; if (typeof Symbol === "function" && typeof Symbol.iterator === "symbol") { _typeof = function _typeof(obj) { return typeof obj; }; } else { _typeof = function _typeof(obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; }; } return _typeof(obj); }

exports.__esModule = true;
exports.default = right;

var C = _interopRequireWildcard(require("./../../../i18n/constants"));

var _utils = require("./../utils");

function _getRequireWildcardCache() { if (typeof WeakMap !== "function") return null; var cache = new WeakMap(); _getRequireWildcardCache = function _getRequireWildcardCache() { return cache; }; return cache; }

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } if (obj === null || _typeof(obj) !== "object" && typeof obj !== "function") { return { default: obj }; } var cache = _getRequireWildcardCache(); if (cache && cache.has(obj)) { return cache.get(obj); } var newObj = {}; var hasPropertyDescriptor = Object.defineProperty && Object.getOwnPropertyDescriptor; for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) { var desc = hasPropertyDescriptor ? Object.getOwnPropertyDescriptor(obj, key) : null; if (desc && (desc.get || desc.set)) { Object.defineProperty(newObj, key, desc); } else { newObj[key] = obj[key]; } } } newObj.default = obj; if (cache) { cache.set(obj, newObj); } return newObj; }

function right(customBordersPlugin) {
  return {
    key: 'borders:right',
    name: function name() {
      var label = this.getTranslatedPhrase(C.CONTEXTMENU_ITEMS_BORDERS_RIGHT);
      var hasBorder = (0, _utils.checkSelectionBorders)(this, 'right');

      if (hasBorder) {
        label = (0, _utils.markSelected)(label);
      }

      return label;
    },
    callback: function callback(key, selected) {
      var hasBorder = (0, _utils.checkSelectionBorders)(this, 'right');
      customBordersPlugin.prepareBorder(selected, 'right', hasBorder);
    }
  };
}