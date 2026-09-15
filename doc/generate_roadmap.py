#!/usr/bin/env python3
"""Generate source/roadmap.rst from two sources:

* Shipped: entries tagged ".. roadmap::" in source/changelog.rst. The
  CHANGELOG is the single source of truth for what shipped and when - to
  make a bullet appear here, add a ".. roadmap::" comment directly beneath
  it (see changelog.rst for existing examples).
* Planned: open GitHub issues labelled "roadmap", fetched from the public
  GitHub REST API (no auth needed, no ``gh`` CLI required).

Nothing else needs to change - re-run this script and rebuild the docs.

Usage:
    python generate_roadmap.py
"""

import json
import re
import urllib.error
import urllib.request
from pathlib import Path

SOURCE = Path(__file__).parent / "source" / "changelog.rst"
OUTPUT = Path(__file__).parent / "source" / "roadmap.rst"

GITHUB_REPO = "openbiosim/biosimspace"
GITHUB_LABEL = "roadmap"

VERSION_LINKED_RE = re.compile(r"^`(?P<version>[\w.]+) <(?P<url>[^>]+)>`_ - (?P<date>.+)$")
VERSION_PLAIN_RE = re.compile(r"^(?P<version>\d[\w.]*) - (?P<date>.+)$")
UNDERLINE_RE = re.compile(r"^-{5,}\s*$")

ATTRIBUTION_RE = re.compile(r"\(`@[\w-]+[^)]*\)")
PR_RE = re.compile(r"\(`#(?P<num>\d+) <(?P<url>[^>]+)>`__\)\.?\s*$")
ROLE_RE = re.compile(r":\w+:`([^<`]+)(?: <[^>]+>)?`")
LITERAL_RE = re.compile(r"``([^`]+)``")
LINK_RE = re.compile(r"`([^`<]+) <([^>]+)>`_+")
BARE_INTERPRETED_RE = re.compile(r"`([^`<>]+)`")


def clean_text(text):
    """Strip RST markup from a changelog bullet and pull out its PR link."""
    pr = None
    match = PR_RE.search(text)
    if match:
        pr = (match.group("num"), match.group("url"))
        text = text[: match.start()].rstrip()

    text = ATTRIBUTION_RE.sub("", text)
    text = ROLE_RE.sub(r"\1", text)
    text = LITERAL_RE.sub(r"<code>\1</code>", text)
    text = LINK_RE.sub(lambda m: f'<a href="{m.group(2)}">{m.group(1)}</a>', text)
    text = BARE_INTERPRETED_RE.sub(r"\1", text)
    text = re.sub(r"\s+", " ", text).strip()
    text = text.rstrip(" .")
    return text, pr


def parse_changelog(path):
    """Return a list of {version, url, date, items} dicts, newest first."""
    lines = path.read_text().splitlines()
    releases = []
    current = None
    i, n = 0, len(lines)

    while i < n:
        line = lines[i]
        m = VERSION_LINKED_RE.match(line)
        m_plain = None if m else VERSION_PLAIN_RE.match(line)

        if (m or m_plain) and i + 1 < n and UNDERLINE_RE.match(lines[i + 1]):
            group = m or m_plain
            current = {
                "version": group.group("version"),
                "url": group.group("url") if m else None,
                "date": group.group("date"),
                "items": [],
            }
            releases.append(current)
            i += 2
            continue

        if current is not None and line.startswith("* "):
            bullet_lines = [line[2:]]
            j = i + 1
            while j < n and lines[j].strip() and not lines[j].startswith(("* ", "..")) \
                    and not (UNDERLINE_RE.match(lines[j])):
                bullet_lines.append(lines[j].strip())
                j += 1
            bullet_text = " ".join(bullet_lines).strip()

            k = j
            while k < n and lines[k].strip() == "":
                k += 1
            tagged = k < n and lines[k].strip() == ".. roadmap::"
            if tagged:
                text, pr = clean_text(bullet_text)
                current["items"].append({"text": text, "pr": pr})
                k += 1
            i = k
            continue

        i += 1

    return [r for r in releases if r["items"]]


def fetch_planned_issues():
    """Fetch open issues labelled GITHUB_LABEL via the public GitHub API.

    Returns an empty list (with a warning) rather than raising, so a
    network hiccup or rate limit never blocks a docs build.
    """
    url = (
        f"https://api.github.com/repos/{GITHUB_REPO}/issues"
        f"?labels={GITHUB_LABEL}&state=open&per_page=100"
    )
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "biosimspace-roadmap-generator",
            "Accept": "application/vnd.github+json",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=10) as response:
            data = json.load(response)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
        print(f"Warning: could not fetch planned roadmap issues from GitHub ({exc}); leaving section empty.")
        return []

    # The issues endpoint also returns pull requests - exclude those.
    return [
        {"number": item["number"], "title": item["title"], "url": item["html_url"]}
        for item in data
        if "pull_request" not in item
    ]


def render_shipped_html(releases):
    rows = []
    for release in releases:
        if release["url"]:
            heading = f'<a href="{release["url"]}">{release["version"]}</a>'
        else:
            heading = release["version"]

        items_html = []
        for item in release["items"]:
            pr_link = ""
            if item["pr"]:
                num, url = item["pr"]
                pr_link = f' <a class="roadmap-pr" href="{url}">#{num}</a>'
            items_html.append(f'<li>{item["text"]}{pr_link}</li>')

        rows.append(f"""
    <div class="roadmap-release">
      <div class="roadmap-marker" aria-hidden="true"></div>
      <div class="roadmap-release-body">
        <div class="roadmap-release-header">
          <span class="roadmap-version">{heading}</span>
          <span class="roadmap-date">{release["date"]}</span>
        </div>
        <ul class="roadmap-items">
          {"".join(items_html)}
        </ul>
      </div>
    </div>""")

    return f"""<div class="roadmap-timeline">
{"".join(rows)}
</div>
"""


def render_planned_html(issues):
    issues_url = (
        f"https://github.com/{GITHUB_REPO}/issues"
        f"?q=is%3Aissue+is%3Aopen+label%3A{GITHUB_LABEL}"
    )
    if not issues:
        return f"""<div class="roadmap-planned">
  <p class="roadmap-planned-empty">
    Nothing tagged yet - see the
    <a href="{issues_url}">open issues labelled "{GITHUB_LABEL}"</a>.
  </p>
</div>
"""

    items_html = "".join(
        f'<li><a href="{issue["url"]}">{issue["title"]}</a> '
        f'<a class="roadmap-pr" href="{issue["url"]}">#{issue["number"]}</a></li>'
        for issue in issues
    )
    return f"""<div class="roadmap-planned">
  <ul class="roadmap-planned-items">
    {items_html}
  </ul>
  <p class="roadmap-planned-footer">
    Tracked via issues labelled
    <a href="{issues_url}">"{GITHUB_LABEL}"</a>.
  </p>
</div>
"""


HEADER = """.. THIS FILE IS AUTO-GENERATED - DO NOT EDIT BY HAND.
.. Run ``python generate_roadmap.py`` from the ``doc`` directory to
.. regenerate it from changelog.rst and open GitHub issues.

Roadmap
=======

A sparse, curated timeline of the capabilities added to :mod:`BioSimSpace`
over time, plus what's currently planned. This isn't every change - see the
:doc:`changelog` for the full, linear history of every fix and tweak - just
the entries that mark a new piece of scientific or engineering capability.

Planned
-------

"""

SHIPPED_HEADER = """
Shipped
-------

"""

STYLE = """
.. raw:: html

    <style>
      .roadmap-planned-items {
        list-style: none;
        margin: 0;
        padding: 0;
      }
      .roadmap-planned-items li {
        padding: 0.5em 0.8em;
        margin: 0.3em 0;
        border: 1px solid var(--color-background-border);
        border-radius: 0.3em;
        background: var(--color-background-secondary);
      }
      .roadmap-planned-empty {
        color: var(--color-foreground-muted);
        font-style: italic;
      }
      .roadmap-planned-footer {
        margin-top: 0.6em;
        font-size: 0.85em;
        color: var(--color-foreground-muted);
      }
      .roadmap-timeline {
        position: relative;
        margin: 2em 0;
        padding-left: 1.6em;
        border-left: 2px solid var(--color-background-border);
      }
      .roadmap-release {
        position: relative;
        margin-bottom: 1.8em;
      }
      .roadmap-marker {
        position: absolute;
        left: -1.95em;
        top: 0.3em;
        width: 0.6em;
        height: 0.6em;
        border-radius: 50%;
        background: var(--color-brand-content);
      }
      .roadmap-release-header {
        display: flex;
        align-items: baseline;
        gap: 0.75em;
        flex-wrap: wrap;
      }
      .roadmap-version {
        font-weight: 600;
        font-size: 1.05em;
      }
      .roadmap-version a {
        color: var(--color-foreground-primary);
        text-decoration: none;
      }
      .roadmap-version a:hover {
        color: var(--color-brand-content);
      }
      .roadmap-date {
        color: var(--color-foreground-muted);
        font-size: 0.85em;
      }
      .roadmap-items {
        margin: 0.4em 0 0 0;
        padding-left: 1.2em;
      }
      .roadmap-items li {
        margin: 0.2em 0;
      }
      .roadmap-items code {
        background: var(--color-background-secondary);
        border-radius: 0.2em;
        padding: 0.05em 0.35em;
      }
      .roadmap-pr {
        font-size: 0.8em;
        color: var(--color-foreground-muted);
        text-decoration: none;
      }
      .roadmap-pr:hover {
        color: var(--color-brand-content);
      }
    </style>
"""


def indent(text, prefix="    "):
    return "\n".join(prefix + line if line else line for line in text.splitlines())


def as_raw_block(html):
    return ".. raw:: html\n\n" + indent(html) + "\n"


def main():
    releases = parse_changelog(SOURCE)
    planned = fetch_planned_issues()

    parts = [
        HEADER,
        as_raw_block(render_planned_html(planned)),
        SHIPPED_HEADER,
        as_raw_block(render_shipped_html(releases)),
        STYLE,
    ]
    OUTPUT.write_text("\n".join(parts))
    print(
        f"Wrote {len(planned)} planned issue(s), "
        f"{len(releases)} release(s) ({sum(len(r['items']) for r in releases)} shipped entries) "
        f"to {OUTPUT}"
    )


if __name__ == "__main__":
    main()
