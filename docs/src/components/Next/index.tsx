import React from "react";
import Link from "@docusaurus/Link";
import styles from "./styles.module.css";

export interface NextItem {
    to: string;
    title: string;
    body?: string;
}

/** A short list of where to go next, shown at the end of a page. */
export default function Next({ items }: { items: NextItem[] }) {
    return (
        <nav className={styles.next} aria-label="Next steps">
            {items.map((item) => (
                <Link key={item.to} to={item.to} className={styles.item}>
                    <span className={styles.title}>{item.title}</span>
                    {item.body && <span className={styles.body}>{item.body}</span>}
                </Link>
            ))}
        </nav>
    );
}
