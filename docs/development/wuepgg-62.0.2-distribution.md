# wuepgg 62.0.2 distribution notice

HelixWeave's genome-browser view pins npm package `wuepgg@62.0.2`.

- Registry tarball: `https://registry.npmjs.org/wuepgg/-/wuepgg-62.0.2.tgz`
- Integrity: `sha512-gXBwD7I5H6YNBXXJJgbcPc1UFglrZ9weGpwU6U1RBnizALUiJMxH3xdUHSZtY/cfmlSQDHhXCXVTOStfx4eOsQ==`
- Package license field: `SEE LICENSE IN LICENSE.txt`
- Upstream text compared on 2026-11-12: `https://raw.githubusercontent.com/twlab/eg3/main/LICENSE.txt`
- The installed `node_modules/wuepgg/LICENSE.txt` matches that upstream text after line-ending and trailing-space normalization. `frontend/src/features/run-artifacts/wuepgg-LICENSE.txt` is the same text. A copy inside the source tree does not by itself put the notice into the built frontend.

The grant excludes commercial use and sales. Copies or substantial portions of the Software must include the agreement and `Copyright (c) 2018 Washington University in St. Louis`. This is not an MIT, Apache, BSD, or other permissive license.

The built frontend carries the agreement at `third-party-notices/wuepgg-62.0.2-LICENSE.html`. Vite copies `frontend/public/` into `dist/`, and the frontend asset packager includes that file in the packaged static set. The packaged app serves `.html` as `text/html`.

`wuepgg@62.0.2` depends on `wuepgg3-track@9.0.3` (`https://registry.npmjs.org/wuepgg3-track/-/wuepgg3-track-9.0.3.tgz`, integrity `sha512-KmyIchdYLLCw4/cZ78PIy0yeJorgsoMBr1Tjz5KmZ772y0QxbBD+dCuoQXVopXlKppB6BDJteDV+p2dcA24BQw==`). The published track package has no `license` field and no LICENSE file, and its distributed JavaScript does not carry a license header. It is not an MIT or Apache grant. The WashU agreement above is the license of the pinned `wuepgg` package.

Other direct `wuepgg` dependencies declare their own SPDX licenses (mostly MIT, plus Apache-2.0 for `firebase` and `fuse.js`, ISC for `qrcode.react`, and `(MIT OR GPL-3.0-or-later)` for `jszip`). Those declarations are not a substitute for the WashU agreement.

HelixWeave's repository license is unchanged. Whether a particular HelixWeave use or distribution is non-commercial, and whether any commercial authorization is required, is an open product decision. This note does not accept that use, contact Washington University, or purchase a license.
