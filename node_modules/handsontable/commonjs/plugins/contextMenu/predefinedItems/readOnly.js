"use strict";

function _typeof(obj) { "@babel/helpers - typeof"; if (typeof Symbol === "function" && typeof Symbol.iterator === "symbol") { _typeof = function _typeof(obj) { return typeof obj; }; } else { _typeof = function _typeof(obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; }; } return _typeof(obj); }

exports.__esModule = true;
exports.default = readOnlyItem;
exports.KEY = void 0;

var _utils = require("./../utils");

var _array = require("./../../../helpers/array");

var C = _interopRequireWildcard(require("./../../../i18n/constants"));

function _getRequireWildcardCache() { if (typeof WeakMap !== "function") return null; var cache = new WeakMap(); _getRequireWildcardCache = function _getRequireWildcardCache() { return cache; }; return cache; }

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } if (obj === null || _typeof(obj) !== "object" && typeof obj !== "function") { return { default: obj }; } var cache = _getRequireWildcardCache(); if (cache && cache.has(obj)) { return cache.get(obj); } var newObj = {}; var hasPropertyDescriptor = Object.defineProperty && Object.getOwnPropertyDescriptor; for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) { var desc = hasPropertyDescriptor ? Object.getOwnPropertyDescriptor(obj, key) : null; if (desc && (desc.get || desc.set)) { Object.defineProperty(newObj, key, desc); } else { newObj[key] = obj[key]; } } } newObj.default = obj; if (cache) { cache.set(obj, newObj); } return newObj; }

var KEY = 'make_read_only';
exports.KEY = KEY;

function readOnlyItem() {
  return {
    key: KEY,
    name: function name() {
      var _this = this;

      var label = this.getTranslatedPhrase(C.CONTEXTMENU_ITEMS_READ_ONLY);
      var atLeastOneReadOnly = (0, _utils.checkSelectionConsistency)(this.getSelectedRange(), function (row, col) {
        return _this.getCellMeta(row, col).readOnly;
      });

      if (atLeastOneReadOnly) {
        label = (0, _utils.markLabelAsSelected)(label);
      }

      return label;
    },
    callback: function callback() {
      var _this2 = this;

      var ranges = this.getSelectedRange();
      var atLeastOneReadOnly = (0, _utils.checkSelectionConsistency)(ranges, function (row, col) {
        return _this2.getCellMeta(row, col).readOnly;
      });
      (0, _array.arrayEach)(ranges, function (range) {
        range.forAll(function (row, col) {
          _this2.setCellMeta(row, col, 'readOnly', !atLeastOneReadOnly);
        });
      });
      this.render();
    },
    disabled: function disabled() {
      return !(this.getSelectedRange() && !this.selection.isSelectedByCorner());
    }
  };
}