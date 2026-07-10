""" allsky_cloudforecast.py

Clear-sky / cloud nowcast module for Allsky.
https://github.com/AllskyTeam/allsky

Author:      Benjamin Hartwich (https://astronomy.garden)
Home / docs: https://github.com/benhartwich/allsky-cloudforecast

Estimates cloud cover from the all-sky image DAY and NIGHT, and turns the recent
trend into a short-term nowcast ("clearing", "clouding over", "stable").

Two complementary measurements, chosen automatically by the Allsky day/night event:

  * Day   - Red/Blue Ratio (RBR): clear sky is blue (low R/B), cloud is white/grey
            (R~=B, high R/B). The standard method used by sky-imager meteorology.
  * Night - star deficit: a clear sky is dotted with stars everywhere; cloud blanks
            them out (template-matched star count on a coarse grid).

A rolling cloud.json feeds a dashboard card. The nowcast is a persistence/trend
extrapolation of the last ~45 min — good for ~30-60 min, not a weather forecast.

Uses only cv2 + numpy (already in the Allsky venv).
"""
import allsky_shared as s
import os
import json
import time
import subprocess
import cv2
import numpy as np

metaData = {
    "name": "Cloud Forecast",
    "description": "Day+night cloud cover from the image (RBR / star deficit) with a short-term clear-sky nowcast",
    "version": "v0.1.0",
    "events": [
        "day",
        "night"
    ],
    "experimental": "false",
    "module": "allsky_cloudforecast",
    "arguments": {
        "mask": "meteor_mask.png",
        "rbr_thr": "0.88",
        "clear_pct": "20",
        "overcast_pct": "65",
        "trend_min": "45",
        "history_hours": "48",
        "publish_web": "true",
        "debug": "false"
    },
    "argumentdetails": {
        "mask": {
            "required": "false",
            "description": "Sky Mask",
            "help": "Mask image (overlay images folder). White = sky to analyse, black = ignore (trees/horizon).",
            "type": {"fieldtype": "image"}
        },
        "rbr_thr": {
            "required": "false",
            "description": "Daytime R/B Cloud Threshold",
            "help": "A sky pixel counts as cloud (daytime) when Red/Blue exceeds this. Clear sky ~0.5-0.6, cloud >0.85.",
            "type": {"fieldtype": "spinner", "min": 0.6, "max": 1.3, "step": 0.01}
        },
        "clear_pct": {
            "required": "false",
            "description": "Clear Below (%)",
            "help": "Cloud cover below this counts as a clear sky",
            "type": {"fieldtype": "spinner", "min": 2, "max": 50, "step": 1}
        },
        "overcast_pct": {
            "required": "false",
            "description": "Overcast Above (%)",
            "help": "Cloud cover above this counts as overcast",
            "type": {"fieldtype": "spinner", "min": 40, "max": 95, "step": 1}
        },
        "trend_min": {
            "required": "false",
            "description": "Nowcast Window (min)",
            "help": "How many minutes of recent history the trend/nowcast is fitted over",
            "type": {"fieldtype": "spinner", "min": 15, "max": 180, "step": 5}
        },
        "history_hours": {
            "required": "false",
            "description": "History (hours)",
            "help": "How much history to keep in cloud.json for charting",
            "type": {"fieldtype": "spinner", "min": 1, "max": 240, "step": 1}
        },
        "publish_web": {
            "required": "false",
            "description": "Publish to Website",
            "help": "Copy cloud.json into the website folder (and upload it) for the dashboard",
            "type": {"fieldtype": "checkbox"}
        },
        "debug": {
            "required": "false",
            "description": "Enable debug images",
            "help": "Write the cloud-mask image to the allsky tmp debug folder",
            "tab": "Debug",
            "type": {"fieldtype": "checkbox"}
        }
    },
    "enabled": "false",
    "changelog": {
        "v0.1.0": [
            {
                "author": "Benjamin Hartwich",
                "authorurl": "https://github.com/benhartwich",
                "changes": "Initial day (RBR) + night (star deficit) cloud cover with a trend-based clear-sky nowcast and dashboard json"
            }
        ]
    }
}

_maskCache = {"name": None, "mask": None}
_starTemplate = None


def _mask(name, shape):
    name = (name or "").strip()
    if not name:
        return np.full(shape, 255, np.uint8)
    if _maskCache["name"] == name and _maskCache["mask"] is not None \
            and _maskCache["mask"].shape == shape:
        return _maskCache["mask"]
    p = os.path.join(s.ALLSKY_OVERLAY, "images", name)
    m = cv2.imread(p, cv2.IMREAD_GRAYSCALE)
    if m is None:
        m = np.full(shape, 255, np.uint8)
    elif m.shape != shape:
        m = cv2.resize(m, (shape[1], shape[0]), interpolation=cv2.INTER_NEAREST)
    _maskCache.update(name=name, mask=m)
    return m


def _cloudDay(bgr, mask, rbr_thr):
    """Daytime cloud fraction via Red/Blue ratio. Overexposed (sun) pixels excluded."""
    b, g, r = cv2.split(bgr.astype(np.float32))
    rbr = r / (b + 1.0)
    sat = (r > 250) & (g > 250) & (b > 250)          # sun / blown highlights
    sky = (mask > 127) & ~sat
    n = int(sky.sum())
    if n < 1000:
        return None, None
    cloud = (rbr > rbr_thr) & sky
    frac = round(100.0 * int(cloud.sum()) / n, 1)
    return frac, (cloud.astype(np.uint8) * 255)


def _starPoints(gray, mask, thr=0.55):
    global _starTemplate
    if _starTemplate is None:
        t = np.zeros((15, 15), np.uint8)
        cv2.circle(t, (7, 7), 3, 255, cv2.FILLED)
        _starTemplate = cv2.blur(t, (2, 2))
    img = cv2.bitwise_and(gray, gray, mask=mask)
    try:
        res = cv2.matchTemplate(img, _starTemplate, cv2.TM_CCOEFF_NORMED)
    except Exception:
        return []
    ys, xs = np.where(res >= thr)
    return list(zip(xs.tolist(), ys.tolist()))


def _cloudNight(gray, mask, cell=80):
    """Night cloud fraction = share of the sky grid with no detected stars."""
    pts = _starPoints(gray, mask)
    h, w = mask.shape
    gw, gh = max(1, w // cell), max(1, h // cell)
    maskC = cv2.resize(mask, (gw, gh), interpolation=cv2.INTER_AREA)
    sky_cells = maskC > 127
    total = int(sky_cells.sum())
    if total == 0:
        return None
    star_grid = np.zeros((gh, gw), bool)
    for x, y in pts:
        star_grid[min(gh - 1, y * gh // h), min(gw - 1, x * gw // w)] = True
    clear = int((sky_cells & star_grid).sum())
    return round(100.0 * (1.0 - clear / total), 1)


def _state(cloud, clear_pct, overcast_pct):
    if cloud is None:
        return "unknown"
    if cloud < clear_pct:
        return "clear"
    if cloud > overcast_pct:
        return "overcast"
    return "partly cloudy"


def _nowcast(history, now_t, now_cloud, trend_min, clear_pct, overcast_pct, method):
    """Linear trend over the last trend_min minutes -> (trend, predicted_30, text).
    Only same-method points are used so the day<->night method switch does not create
    a spurious jump, and a minimum time span is required to avoid a degenerate slope."""
    cutoff = now_t - trend_min * 60
    pts = [(d["t"], d["cloud"]) for d in history
           if d.get("t", 0) >= cutoff and d.get("cloud") is not None
           and d.get("method") == method]
    pts.append((now_t, now_cloud))
    span_min = (now_t - min(t for t, _ in pts)) / 60.0
    if len(pts) < 3 or span_min < max(10.0, trend_min * 0.3):
        return "unknown", None, "building history"
    tm = np.array([(t - now_t) / 60.0 for t, _ in pts])     # minutes, <=0
    cl = np.array([c for _, c in pts], float)
    slope = float(np.polyfit(tm, cl, 1)[0])                  # %/min
    pred = float(np.clip(now_cloud + slope * 30.0, 0, 100))
    if slope < -0.4:
        trend = "clearing"
    elif slope > 0.4:
        trend = "clouding over"
    else:
        trend = "stable"
    # short verbal nowcast
    if trend == "clearing" and now_cloud >= clear_pct:
        mins = (now_cloud - clear_pct) / max(1e-3, -slope)
        text = f"clearing — likely clear within ~{int(round(mins/5)*5)} min" if mins < 180 else "slowly clearing"
    elif trend == "clouding over" and now_cloud <= overcast_pct:
        mins = (overcast_pct - now_cloud) / max(1e-3, slope)
        text = f"clouding over — overcast in ~{int(round(mins/5)*5)} min" if mins < 180 else "slowly clouding over"
    else:
        text = f"{_state(now_cloud, clear_pct, overcast_pct)}, {trend}"
    return trend, round(pred, 1), text


def _websiteDir():
    website = s.getEnvironmentVariable("ALLSKY_WEBSITE")
    if not website:
        website = os.path.join(s.getEnvironmentVariable("ALLSKY_HOME") or os.path.expanduser("~/allsky"),
                               "html", "allsky")
    return website


def _uploadRemote(local, fname):
    try:
        if s.getSetting("useremotewebsite") != "true":
            return
        scripts = s.getEnvironmentVariable("ALLSKY_SCRIPTS") or \
            os.path.join(s.getEnvironmentVariable("ALLSKY_HOME") or os.path.expanduser("~/allsky"), "scripts")
        uploader = os.path.join(scripts, "upload.sh")
        if not os.path.isfile(uploader) or not os.path.isfile(local):
            return
        rdir = (s.getSetting("remotewebsiteimagedir") or "").rstrip("/")
        subprocess.Popen([uploader, "--silent", "--wait", "--remote-web", local, rdir, fname, "CloudForecast"],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as ex:
        s.log(1, f"WARNING: cloudforecast remote upload failed: {ex}")


def _appendHistory(record, hours, publish_web):
    path = os.path.join(s.ALLSKY_TMP, "cloud.json")
    try:
        data = json.load(open(path)) if os.path.exists(path) else []
    except Exception:
        data = []
    data.append(record)
    cutoff = record["t"] - hours * 3600
    data = [d for d in data if d.get("t", 0) >= cutoff][-6000:]
    try:
        json.dump(data, open(path, "w"))
    except Exception as ex:
        s.log(1, f"WARNING: cloudforecast could not write history: {ex}")
        return data
    if publish_web:
        try:
            ddir = _websiteDir()
            os.makedirs(ddir, exist_ok=True)
            webpath = os.path.join(ddir, "cloud.json")
            json.dump(data, open(webpath, "w"))
            _uploadRemote(webpath, "cloud.json")
        except Exception as ex:
            s.log(1, f"WARNING: cloudforecast could not publish: {ex}")
    return data


def cloudforecast(params, event):
    if s.image is None:
        return "No image available"

    rbr_thr = s.asfloat(params.get("rbr_thr", 0.88))
    clear_pct = s.asfloat(params.get("clear_pct", 20))
    overcast_pct = s.asfloat(params.get("overcast_pct", 65))
    trend_min = s.int(params.get("trend_min", 45))
    debug = params.get("debug", False)

    shape = s.image.shape[:2]
    mask = _mask(params.get("mask", ""), shape)

    if event == "day":
        cloud, cmask = _cloudDay(s.image, mask, rbr_thr)
        method = "rbr"
    else:
        gray = cv2.cvtColor(s.image, cv2.COLOR_BGR2GRAY) if len(s.image.shape) == 3 else s.image
        cloud = _cloudNight(gray, mask)
        cmask = None
        method = "stars"

    if cloud is None:
        return "Could not measure cloud cover"
    if debug and cmask is not None:
        s.startModuleDebug(metaData["module"])
        s.writeDebugImage(metaData["module"], "cloudmask.png", cmask)

    now = int(time.time())
    history = []
    path = os.path.join(s.ALLSKY_TMP, "cloud.json")
    try:
        history = json.load(open(path)) if os.path.exists(path) else []
    except Exception:
        history = []

    state = _state(cloud, clear_pct, overcast_pct)
    trend, pred30, text = _nowcast(history, now, cloud, trend_min, clear_pct, overcast_pct, method)

    rec = {"t": now, "cloud": cloud, "method": method, "state": state,
           "trend": trend, "pred30": pred30}
    _appendHistory(rec, s.int(params.get("history_hours", 48)), params.get("publish_web", True))

    s.setEnvironmentVariable("AS_CLOUDFRAC", f"{cloud:.0f}")
    s.setEnvironmentVariable("AS_SKYSTATE", state)
    s.setEnvironmentVariable("AS_CLOUDTREND", trend)
    s.setEnvironmentVariable("AS_CLOUDNOWCAST", text)

    result = f"Cloud {cloud:.0f}% ({method}) — {state}; {text}"
    s.log(4, f"INFO: {result}")
    return result


def cloudforecast_cleanup():
    moduleData = {
        "metaData": metaData,
        "cleanup": {
            "files": {os.path.join(s.ALLSKY_TMP, "cloud.json")},
            "env": {"AS_CLOUDFRAC", "AS_SKYSTATE", "AS_CLOUDTREND", "AS_CLOUDNOWCAST"}
        }
    }
    s.cleanupModule(moduleData)
