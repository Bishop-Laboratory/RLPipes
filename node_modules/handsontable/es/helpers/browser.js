import { objectEach } from './object';

var tester = function tester(testerFunc) {
  var result = {
    value: false
  };

  result.test = function (ua, vendor) {
    result.value = testerFunc(ua, vendor);
  };

  return result;
};

var browsers = {
  chrome: tester(function (ua, vendor) {
    return /Chrome/.test(ua) && /Google/.test(vendor);
  }),
  edge: tester(function (ua) {
    return /Edge/.test(ua);
  }),
  firefox: tester(function (ua) {
    return /Firefox/.test(ua);
  }),
  ie: tester(function (ua) {
    return /Trident/.test(ua);
  }),
  // eslint-disable-next-line no-restricted-globals
  ie8: tester(function () {
    return !document.createTextNode('test').textContent;
  }),
  // eslint-disable-next-line no-restricted-globals
  ie9: tester(function () {
    return !!document.documentMode;
  }),
  mobile: tester(function (ua) {
    return /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(ua);
  }),
  safari: tester(function (ua, vendor) {
    return /Safari/.test(ua) && /Apple Computer/.test(vendor);
  })
};
var platforms = {
  mac: tester(function (platform) {
    return /^Mac/.test(platform);
  }),
  win: tester(function (platform) {
    return /^Win/.test(platform);
  }),
  linux: tester(function (platform) {
    return /^Linux/.test(platform);
  })
};
/**
 * @param {object} [metaObject] The browser identity collection.
 * @param {object} [metaObject.userAgent] The user agent reported by browser.
 * @param {object} [metaObject.vendor] The vendor name reported by browser.
 */

export function setBrowserMeta() {
  var _ref = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {},
      _ref$userAgent = _ref.userAgent,
      userAgent = _ref$userAgent === void 0 ? navigator.userAgent : _ref$userAgent,
      _ref$vendor = _ref.vendor,
      vendor = _ref$vendor === void 0 ? navigator.vendor : _ref$vendor;

  objectEach(browsers, function (_ref2) {
    var test = _ref2.test;
    return void test(userAgent, vendor);
  });
}
/**
 * @param {object} [metaObject] The platform identity collection.
 * @param {object} [metaObject.platform] The platform ID.
 */

export function setPlatformMeta() {
  var _ref3 = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {},
      _ref3$platform = _ref3.platform,
      platform = _ref3$platform === void 0 ? navigator.platform : _ref3$platform;

  objectEach(platforms, function (_ref4) {
    var test = _ref4.test;
    return void test(platform);
  });
}
setBrowserMeta();
setPlatformMeta();
export function isChrome() {
  return browsers.chrome.value;
}
export function isEdge() {
  return browsers.edge.value;
}
export function isIE() {
  return browsers.ie.value;
}
export function isIE8() {
  return browsers.ie8.value;
}
export function isIE9() {
  return browsers.ie9.value;
}
export function isMSBrowser() {
  return browsers.ie.value || browsers.edge.value;
}
export function isMobileBrowser() {
  return browsers.mobile.value;
}
export function isSafari() {
  return browsers.safari.value;
}
export function isFirefox() {
  return browsers.firefox.value;
}
/**
 * @returns {boolean}
 */

export function isWindowsOS() {
  return platforms.win.value;
}
/**
 * @returns {boolean}
 */

export function isMacOS() {
  return platforms.mac.value;
}
/**
 * @returns {boolean}
 */

export function isLinuxOS() {
  return platforms.linux.value;
}