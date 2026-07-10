# Changelog

## v0.1.0
- Initial release.
- Daytime cloud cover via Red/Blue Ratio; night-time via star deficit; method chosen
  automatically by the Allsky day/night event.
- Trend-based clear-sky nowcast (clearing / clouding over / stable) with a short-term
  estimate, using same-method points and a minimum-time-span guard.
- Rolling `cloud.json` for charting; overlay environment variables.
