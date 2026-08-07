#!/usr/bin/env python3
"""
Collect OMIM update entries from either:
1) OMIM monthly update web pages, e.g. https://www.omim.org/statistics/updates/2026/1
2) saved OMIM notification emails as .txt/.eml files
3) an IMAP mailbox, filtering messages from do-not-reply@omim.org

Output: a normalized TSV of OMIM update entries.
"""

import argparse
import csv
import datetime as dt
import email
import html
import imaplib
import os
import re
import sys
import urllib.request
from email.header import decode_header
from email.message import Message
from typing import Dict, Iterable, List, Optional

DATE_RE = re.compile(r"^[A-Z][a-z]+\s+\d{1,2}(?:st|nd|rd|th),\s+\d{4}$")
ENTRY_RE = re.compile(r"^([#*+%])?(\d{6})(.+)$")
EMAIL_ENTRY_RE = re.compile(r"^(\d{6})\s+(.+)$")
GENE_HINT_RE = re.compile(r"\b(?:in|of|for)\s+the\s+([A-Z][A-Z0-9-]+)\s+gene\b|\b([A-Z][A-Z0-9-]+)\s+gene\b")

DEFAULT_SECTIONS = ["New Entries", "New Clinical Synopses"]
ALL_SECTIONS = ["New Entries", "New Clinical Synopses", "Updated Clinical Synopses", "Updated Entries"]


def month_range(start_ym: str, end_ym: str):
    y, m = map(int, start_ym.split("-"))
    ey, em = map(int, end_ym.split("-"))
    cur = dt.date(y, m, 1)
    end = dt.date(ey, em, 1)
    while cur <= end:
        yield cur.year, cur.month
        if cur.month == 12:
            cur = dt.date(cur.year + 1, 1, 1)
        else:
            cur = dt.date(cur.year, cur.month + 1, 1)


def default_last_quarter(today: Optional[dt.date] = None):
    if today is None:
        today = dt.date.today()
    # previous completed quarter by default
    q = (today.month - 1) // 3 + 1
    if q == 1:
        year = today.year - 1
        start_m, end_m = 10, 12
    elif q == 2:
        year = today.year
        start_m, end_m = 1, 3
    elif q == 3:
        year = today.year
        start_m, end_m = 4, 6
    else:
        year = today.year
        start_m, end_m = 7, 9
    return f"{year}-{start_m:02d}", f"{year}-{end_m:02d}"


def fetch_url(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Mozilla/5.0 OMIM-reanalysis-pipeline/1.0",
            "Accept": "text/html,text/plain;q=0.9,*/*;q=0.8",
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = r.read()
    return data.decode("utf-8", errors="replace")


def html_to_text(raw: str) -> str:
    # Preserve line breaks for list-like pages before tag removal.
    x = re.sub(r"(?i)<\s*br\s*/?>", "\n", raw)
    x = re.sub(r"(?i)</\s*(p|div|li|h\d|tr|table|section)\s*>", "\n", x)
    x = re.sub(r"<[^>]+>", "", x)
    x = html.unescape(x)
    lines = [re.sub(r"\s+", " ", line).strip() for line in x.splitlines()]
    return "\n".join([line for line in lines if line])


def parse_update_text(text: str, source: str, url: str = "", keep_sections: Optional[List[str]] = None) -> List[Dict[str, str]]:
    if keep_sections is None:
        keep_sections = DEFAULT_SECTIONS
    keep_set = set(s.rstrip(":") for s in keep_sections)
    known_section_set = set(s.rstrip(":") for s in ALL_SECTIONS)

    current_date = ""
    current_section = ""
    rows: List[Dict[str, str]] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if DATE_RE.match(line):
            current_date = line
            current_section = ""
            continue

        if line.endswith(":"):
            sec = line.rstrip(":").strip()
            if sec in known_section_set:
                current_section = sec
            else:
                current_section = ""
            continue

        if current_section not in keep_set:
            continue

        m = ENTRY_RE.match(line)
        if not m:
            continue

        prefix, mim_number, title = m.group(1) or "", m.group(2), m.group(3).strip()
        rows.append({
            "MIM_number": mim_number,
            "Title": title,
            "Date": current_date,
            "Section": current_section,
            "Prefix": prefix,
            "Source": source,
            "URL": url,
            "Email_gene_hint": "",
        })

    return rows


def decode_header_value(x: str) -> str:
    parts = decode_header(x or "")
    out = []
    for val, enc in parts:
        if isinstance(val, bytes):
            out.append(val.decode(enc or "utf-8", errors="replace"))
        else:
            out.append(val)
    return "".join(out)


def message_to_text(msg: Message) -> str:
    parts = []
    if msg.is_multipart():
        for part in msg.walk():
            ctype = part.get_content_type()
            disp = str(part.get("Content-Disposition", ""))
            if "attachment" in disp.lower():
                continue
            if ctype in ("text/plain", "text/html"):
                payload = part.get_payload(decode=True)
                if payload is None:
                    continue
                charset = part.get_content_charset() or "utf-8"
                text = payload.decode(charset, errors="replace")
                if ctype == "text/html":
                    text = html_to_text(text)
                parts.append(text)
    else:
        payload = msg.get_payload(decode=True)
        if payload is not None:
            charset = msg.get_content_charset() or "utf-8"
            text = payload.decode(charset, errors="replace")
            if msg.get_content_type() == "text/html":
                text = html_to_text(text)
            parts.append(text)
    return "\n".join(parts)


def parse_email_text(text: str, source: str = "email", msg_date: str = "") -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    gene_hint = ""
    for g in GENE_HINT_RE.findall(text):
        gene_hint = g[0] or g[1]
        if gene_hint:
            break

    # Email body usually has lines like:
    # 621656 RETINITIS PIGMENTOSA 109; RP109
    for raw_line in text.splitlines():
        line = raw_line.strip()
        m = EMAIL_ENTRY_RE.match(line)
        if not m:
            continue
        mim_number, title = m.group(1), m.group(2).strip()
        # Keep only plausible OMIM entry title lines, not phone/date lines.
        if len(title) < 3 or title.lower().startswith("unsubscribe"):
            continue
        rows.append({
            "MIM_number": mim_number,
            "Title": title,
            "Date": msg_date,
            "Section": "Email New gene/phenotype relationship",
            "Prefix": "#",
            "Source": source,
            "URL": "",
            "Email_gene_hint": gene_hint,
        })
    return rows


def parse_email_file(path: str) -> List[Dict[str, str]]:
    with open(path, "rb") as f:
        data = f.read()
    try:
        msg = email.message_from_bytes(data)
        subject = decode_header_value(msg.get("Subject", ""))
        msg_date = msg.get("Date", "")
        text = message_to_text(msg)
        if not text.strip():
            text = data.decode("utf-8", errors="replace")
        return parse_email_text(text, source=f"email_file:{os.path.basename(path)}", msg_date=msg_date)
    except Exception:
        text = data.decode("utf-8", errors="replace")
        return parse_email_text(text, source=f"email_text:{os.path.basename(path)}")


def collect_from_imap(args) -> List[Dict[str, str]]:
    password = os.environ.get(args.imap_password_env, "")
    if not password:
        raise RuntimeError(f"Environment variable {args.imap_password_env} is empty. Do not put passwords in scripts.")
    imap = imaplib.IMAP4_SSL(args.imap_server, args.imap_port)
    imap.login(args.imap_user, password)
    imap.select(args.imap_folder)

    criteria = ['FROM', '"do-not-reply@omim.org"']
    if args.imap_since:
        criteria += ['SINCE', args.imap_since]
    status, data = imap.search(None, *criteria)
    if status != "OK":
        raise RuntimeError(f"IMAP search failed: {status}")

    rows: List[Dict[str, str]] = []
    ids = data[0].split()
    for msg_id in ids:
        status, msg_data = imap.fetch(msg_id, "(RFC822)")
        if status != "OK":
            continue
        raw = msg_data[0][1]
        msg = email.message_from_bytes(raw)
        subject = decode_header_value(msg.get("Subject", ""))
        if args.subject_contains and args.subject_contains.lower() not in subject.lower():
            continue
        msg_date = msg.get("Date", "")
        text = message_to_text(msg)
        rows.extend(parse_email_text(text, source="imap", msg_date=msg_date))
    imap.logout()
    return rows


def dedup_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    seen = set()
    out = []
    for r in rows:
        key = (r.get("MIM_number", ""), r.get("Date", ""), r.get("Section", ""), r.get("Email_gene_hint", ""))
        if key in seen:
            continue
        seen.add(key)
        out.append(r)
    return out


def main():
    p = argparse.ArgumentParser(description="Collect OMIM updated disease entries for reanalysis.")
    p.add_argument("--start", help="Start month YYYY-MM. Default: start of previous completed quarter.")
    p.add_argument("--end", help="End month YYYY-MM. Default: end of previous completed quarter.")
    p.add_argument("--sections", default=",".join(DEFAULT_SECTIONS), help="Comma-separated sections to keep.")
    p.add_argument("--web", action="store_true", help="Collect OMIM monthly update pages.")
    p.add_argument("--email-files", nargs="*", default=[], help="Saved OMIM notification emails/text files.")
    p.add_argument("--imap", action="store_true", help="Collect OMIM notification emails via IMAP.")
    p.add_argument("--imap-server", default="imap.163.com")
    p.add_argument("--imap-port", type=int, default=993)
    p.add_argument("--imap-user", default="")
    p.add_argument("--imap-folder", default="INBOX")
    p.add_argument("--imap-password-env", default="OMIM_MAIL_PASSWORD")
    p.add_argument("--imap-since", default="", help="IMAP SINCE date, e.g. 01-Jan-2026.")
    p.add_argument("--subject-contains", default="OMIM")
    p.add_argument("--out", required=True)
    p.add_argument("--raw-dir", default="omim_raw_pages")
    args = p.parse_args()

    if not args.start or not args.end:
        s, e = default_last_quarter()
        args.start = args.start or s
        args.end = args.end or e

    keep_sections = [x.strip() for x in args.sections.split(",") if x.strip()]
    rows: List[Dict[str, str]] = []

    if args.web:
        os.makedirs(args.raw_dir, exist_ok=True)
        for year, month in month_range(args.start, args.end):
            url = f"https://www.omim.org/statistics/updates/{year}/{month}"
            print(f"[INFO] fetching {url}", file=sys.stderr)
            try:
                raw = fetch_url(url)
            except Exception as e:
                print(f"[WARN] failed to fetch {url}: {e}", file=sys.stderr)
                continue
            raw_path = os.path.join(args.raw_dir, f"omim_updates_{year}_{month:02d}.html")
            with open(raw_path, "w", encoding="utf-8") as f:
                f.write(raw)
            text = html_to_text(raw)
            with open(raw_path.replace(".html", ".txt"), "w", encoding="utf-8") as f:
                f.write(text)
            rows.extend(parse_update_text(text, source="web", url=url, keep_sections=keep_sections))

    for path in args.email_files:
        rows.extend(parse_email_file(path))

    if args.imap:
        rows.extend(collect_from_imap(args))

    rows = dedup_rows(rows)
    fieldnames = ["MIM_number", "Title", "Date", "Section", "Prefix", "Source", "URL", "Email_gene_hint"]
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    print(f"[DONE] wrote {len(rows)} records to {args.out}")


if __name__ == "__main__":
    main()
