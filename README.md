# allsky_cloudforecast

A **day + night cloud cover and clear-sky nowcast** module for
[Allsky](https://github.com/AllskyTeam/allsky).

It estimates how cloudy the sky is directly from the all-sky image — around the
clock — and turns the recent trend into a short-term nowcast: *clearing*,
*clouding over* or *stable*.

## Two methods, chosen automatically

| Time | Method | Idea |
|---|---|---|
| **Day** | Red/Blue Ratio (RBR) | Clear sky is blue (low R/B); cloud is white/grey (R≈B, high R/B). A sky pixel counts as cloud when R/B exceeds a threshold. This is the standard sky-imager meteorology method. |
| **Night** | Star deficit | A clear sky is dotted with stars everywhere; cloud blanks them out. Cloud cover = share of the sky grid with no detected stars. |

The Allsky day/night event picks the method, so a single number tracks cloud cover
continuously through the day.

## Nowcast

The last ~45 minutes of cloud cover are fitted with a straight line:

- falling → **clearing** (with an estimated "clear within ~N min")
- rising → **clouding over** (with "overcast in ~N min")
- flat → **stable**

Only same-method points are used, so the day↔night switch never creates a spurious
jump, and a minimum time span is required before a trend is reported.

> **Honest scope.** This is a persistence/trend extrapolation — reliable for roughly
> the next 30–60 minutes, not a weather forecast. Cloud-motion tracking (optical flow)
> for directional nowcasting is a future improvement.

## Installation

```bash
cp allsky_cloudforecast.py ~/allsky/scripts/modules/
```

Enable **"Cloud Forecast"** in the Allsky WebUI for **both** the day and night flows.
Uses only `cv2` + `numpy` (already in the Allsky venv).

## Configuration

| Setting | Default | Meaning |
|---|---|---|
| Sky Mask | `meteor_mask.png` | White = sky to analyse, black = ignore (trees/horizon) |
| Daytime R/B Cloud Threshold | `0.88` | R/B above which a daytime pixel is cloud (clear ≈ 0.5–0.6) |
| Clear Below (%) | `20` | Cloud cover below this = clear sky |
| Overcast Above (%) | `65` | Cloud cover above this = overcast |
| Nowcast Window (min) | `45` | History window the trend is fitted over |
| History (hours) | `48` | How much history to keep in `cloud.json` |
| Publish to Website | on | Copy/upload `cloud.json` for the dashboard |

The daytime threshold and the clear/overcast bands are site-dependent (light
pollution, horizon glow); tune them once against a clear and an overcast frame.

## Output

- Environment variables `AS_CLOUDFRAC`, `AS_SKYSTATE`, `AS_CLOUDTREND`,
  `AS_CLOUDNOWCAST` — usable in the Allsky overlay.
- A rolling **`cloud.json`**: one `{t, cloud, method, state, trend, pred30}` record
  per frame, ready for a dashboard chart.

## Credits

- [Allsky](https://github.com/AllskyTeam/allsky) by Thomas Jacquin and team.
- Built for [astronomy.garden](https://astronomy.garden).

## License

MIT — see [LICENSE](LICENSE).
